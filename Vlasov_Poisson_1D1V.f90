
!===============================================================
! Vlasov_Poisson_1D1V.f90
!
! 1D-1V electrostatic Vlasov-Poisson solver
!
! Model:
!   df/dt + d(v*f)/dx + d(a*f)/dv = 0
!   dE/dx = rho/epsilon0
!   rho = q*(n - n0),  n = integral f dv
!   a = (q/m)*E
!
! Normalized units:
!   m = 1, |q| = 1, epsilon0 = 1, n0 = 1
!   q/m = -1 (electron-like species)
!
! Boundary conditions:
!   x: periodic
!   v: zero particle flux at the velocity boundaries
!
! Numerical method:
!   First-order conservative upwind finite-volume fluxes
!   Adaptive timestep based on a combined CFL estimate
!
! Outputs:
!   diagnostics.csv
!   f_final.dat
!===============================================================

module vlasov_poisson_module
    use iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: run_simulation

contains

    subroutine solve_electric_field(f, v, dx, dv, nx, nv, q, n0, eps0, &
                                    density, rho, electric_field)
        implicit none

        integer, intent(in) :: nx, nv
        real(dp), intent(in) :: f(nx,nv), v(nv), dx, dv
        real(dp), intent(in) :: q, n0, eps0
        real(dp), intent(out) :: density(nx), rho(nx), electric_field(nx)

        real(dp) :: eface(0:nx), mean_eface
        integer :: i, j

        ! Number density: n(x) = integral f(x,v) dv
        density = 0.0_dp
        do j = 1, nv
            density(:) = density(:) + f(:,j)*dv
        end do

        ! Uniform immobile background neutralizes the mean charge.
        rho = q*(density - n0)

        ! Remove only the residual mean charge to ensure periodic
        ! compatibility of the discrete Poisson equation.
        rho = rho - sum(rho)/real(nx,dp)

        ! Integrate Gauss's law from cell-centered charge to faces:
        ! E(i+1/2) - E(i-1/2) = rho(i)*dx/eps0
        eface(0) = 0.0_dp
        do i = 1, nx
            eface(i) = eface(i-1) + rho(i)*dx/eps0
        end do

        ! The periodic electric field is defined up to a constant.
        ! Set its spatial mean to zero.
        mean_eface = sum(eface(0:nx-1))/real(nx,dp)
        eface = eface - mean_eface

        ! Cell-centered electric field.
        do i = 1, nx
            electric_field(i) = 0.5_dp*(eface(i-1) + eface(i))
        end do
    end subroutine solve_electric_field


    subroutine calculate_diagnostics(f, v, electric_field, dx, dv, nx, nv, &
                                     m, eps0, mass, kinetic_energy, &
                                     field_energy)
        implicit none

        integer, intent(in) :: nx, nv
        real(dp), intent(in) :: f(nx,nv), v(nv), electric_field(nx)
        real(dp), intent(in) :: dx, dv, m, eps0
        real(dp), intent(out) :: mass, kinetic_energy, field_energy

        integer :: i, j

        mass = 0.0_dp
        kinetic_energy = 0.0_dp

        do j = 1, nv
            mass = mass + sum(f(:,j))*dx*dv
            kinetic_energy = kinetic_energy + &
                0.5_dp*m*v(j)**2*sum(f(:,j))*dx*dv
        end do

        field_energy = 0.5_dp*eps0*sum(electric_field**2)*dx
    end subroutine calculate_diagnostics


    subroutine write_final_distribution(f, x, v, nx, nv)
        implicit none

        integer, intent(in) :: nx, nv
        real(dp), intent(in) :: f(nx,nv), x(nx), v(nv)
        integer :: i, j, unit_number

        open(newunit=unit_number, file="f_final.dat", status="replace", &
             action="write")

        write(unit_number,'(a)') "# x  v  f(x,v)"

        do j = 1, nv
            do i = 1, nx
                write(unit_number,'(3(es24.16,1x))') x(i), v(j), f(i,j)
            end do
            write(unit_number,*)
        end do

        close(unit_number)
    end subroutine write_final_distribution


    subroutine run_simulation()
        implicit none

        !---------------- Grid and simulation settings ----------------
        integer, parameter :: nx = 128
        integer, parameter :: nv = 128
        integer, parameter :: max_steps = 200000
        integer, parameter :: diagnostic_interval = 100

        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), parameter :: length_x = 2.0_dp*pi
        real(dp), parameter :: vmin = -6.0_dp
        real(dp), parameter :: vmax =  6.0_dp
        real(dp), parameter :: n0 = 1.0_dp
        real(dp), parameter :: q = -1.0_dp
        real(dp), parameter :: m = 1.0_dp
        real(dp), parameter :: eps0 = 1.0_dp
        real(dp), parameter :: perturbation = 0.01_dp
        real(dp), parameter :: cfl = 0.45_dp
        real(dp), parameter :: t_end = 10.0_dp

        !---------------- Arrays and variables ------------------------
        real(dp) :: dx, dv, time, dt, acceleration_max
        real(dp) :: mass, kinetic_energy, field_energy, total_energy
        real(dp) :: initial_mass, initial_energy
        real(dp) :: velocity_norm, wave_number
        real(dp) :: face_velocity, acceleration
        real(dp) :: flux_left, flux_right
        real(dp) :: mass_error, energy_error
        real(dp) :: max_field, max_acceleration
        real(dp), allocatable :: x(:), v(:)
        real(dp), allocatable :: f(:,:), fnew(:,:)
        real(dp), allocatable :: fx(:,:), fv(:,:)
        real(dp), allocatable :: density(:), rho(:), electric_field(:)

        integer :: i, j, ip, im, step, unit_diag
        integer :: diagnostic_status

        allocate(x(nx), v(nv))
        allocate(f(nx,nv), fnew(nx,nv))
        allocate(fx(nx,nv), fv(nx,0:nv))
        allocate(density(nx), rho(nx), electric_field(nx))

        dx = length_x/real(nx,dp)
        dv = (vmax-vmin)/real(nv,dp)
        wave_number = 2.0_dp*pi/length_x

        do i = 1, nx
            x(i) = (real(i,dp)-0.5_dp)*dx
        end do

        do j = 1, nv
            v(j) = vmin + (real(j,dp)-0.5_dp)*dv
        end do

        !---------------- Initial distribution -----------------------
        ! Maxwellian in velocity, with a small sinusoidal density
        ! perturbation in space.
        !
        ! Normalize the discrete velocity Maxwellian so that its
        ! numerical integral equals one.

        velocity_norm = sum(exp(-0.5_dp*v**2))*dv

        do i = 1, nx
            do j = 1, nv
                f(i,j) = n0 * &
                    (1.0_dp + perturbation*cos(wave_number*x(i))) * &
                    exp(-0.5_dp*v(j)**2) / velocity_norm
            end do
        end do

        call solve_electric_field(f, v, dx, dv, nx, nv, q, n0, eps0, &
                                  density, rho, electric_field)

        call calculate_diagnostics(f, v, electric_field, dx, dv, nx, nv, &
                                   m, eps0, mass, kinetic_energy, field_energy)

        initial_mass = mass
        initial_energy = kinetic_energy + field_energy

        !---------------- Diagnostics file ---------------------------
        open(newunit=unit_diag, file="diagnostics.csv", status="replace", &
             action="write")

        write(unit_diag,'(a)') &
            "step,time,mass,kinetic_energy,electric_field_energy," // &
            "total_energy,mass_relative_error,energy_relative_error,max_abs_E"

        write(unit_diag,'(i0,8(",",es20.12))') 0, 0.0_dp, mass, &
            kinetic_energy, field_energy, kinetic_energy+field_energy, &
            0.0_dp, 0.0_dp, maxval(abs(electric_field))

        write(*,'(a)') "Starting 1D-1V Vlasov-Poisson simulation"
        write(*,'(a,i0,a,i0)') "Grid: nx = ", nx, ", nv = ", nv
        write(*,'(a,es12.4)') "Initial mass: ", initial_mass
        write(*,'(a,es12.4)') "Initial total energy: ", initial_energy

        !---------------- Main time loop -----------------------------
        time = 0.0_dp
        step = 0

        do while (time < t_end .and. step < max_steps)

            ! Recompute E from the current distribution.
            call solve_electric_field(f, v, dx, dv, nx, nv, q, n0, eps0, &
                                      density, rho, electric_field)

            max_field = maxval(abs(electric_field))
            max_acceleration = abs(q/m)*max_field

            ! Combined CFL estimate for x- and v-advection.
            acceleration_max = maxval(abs(v))/dx + max_acceleration/dv

            if (acceleration_max <= tiny(1.0_dp)) then
                dt = t_end-time
            else
                dt = cfl/acceleration_max
            end if

            dt = min(dt, t_end-time)

            !---------------- Spatial fluxes -------------------------
            ! fx(i,j) is the flux through the right face of cell i.
            ! Periodic boundary: cell nx is followed by cell 1.
            do i = 1, nx
                ip = mod(i,nx) + 1

                do j = 1, nv
                    face_velocity = v(j)

                    if (face_velocity >= 0.0_dp) then
                        fx(i,j) = face_velocity*f(i,j)
                    else
                        fx(i,j) = face_velocity*f(ip,j)
                    end if
                end do
            end do

            !---------------- Velocity-space fluxes -----------------
            ! Zero flux at the two velocity boundaries.
            fv(:,0) = 0.0_dp
            fv(:,nv) = 0.0_dp

            do i = 1, nx
                acceleration = (q/m)*electric_field(i)

                do j = 1, nv-1
                    if (acceleration >= 0.0_dp) then
                        fv(i,j) = acceleration*f(i,j)
                    else
                        fv(i,j) = acceleration*f(i,j+1)
                    end if
                end do
            end do

            !---------------- Conservative update -------------------
            do i = 1, nx
                im = i - 1
                if (im < 1) im = nx

                do j = 1, nv
                    flux_left = fx(im,j)
                    flux_right = fx(i,j)

                    fnew(i,j) = f(i,j) &
                        - dt/dx*(flux_right-flux_left) &
                        - dt/dv*(fv(i,j)-fv(i,j-1))
                end do
            end do

            f = fnew
            time = time + dt
            step = step + 1

            !---------------- Diagnostics ----------------------------
            if (mod(step,diagnostic_interval) == 0 .or. time >= t_end) then

                call solve_electric_field(f, v, dx, dv, nx, nv, q, n0, eps0, &
                                          density, rho, electric_field)

                call calculate_diagnostics(f, v, electric_field, dx, dv, &
                                           nx, nv, m, eps0, mass, &
                                           kinetic_energy, field_energy)

                total_energy = kinetic_energy + field_energy

                mass_error = abs(mass-initial_mass) / &
                             max(abs(initial_mass), tiny(1.0_dp))

                energy_error = abs(total_energy-initial_energy) / &
                               max(abs(initial_energy), tiny(1.0_dp))

                write(unit_diag,'(i0,8(",",es20.12))') step, time, mass, &
                    kinetic_energy, field_energy, total_energy, &
                    mass_error, energy_error, maxval(abs(electric_field))

                write(*,'(a,i8,a,f9.4,a,es11.3,a,es11.3)') &
                    "Step ", step, "  t=", time, &
                    "  mass error=", mass_error, &
                    "  energy error=", energy_error
            end if

            ! Stop early if the distribution becomes invalid.
            if (.not. all(f == f)) then
                write(*,'(a)') "ERROR: NaN detected in distribution."
                exit
            end if

            if (maxval(abs(f)) > 1.0e100_dp) then
                write(*,'(a)') "ERROR: distribution magnitude is unbounded."
                exit
            end if

        end do

        close(unit_diag)

        call write_final_distribution(f, x, v, nx, nv)

        write(*,'(a,i0)') "Finished at step: ", step
        write(*,'(a,es12.4)') "Final simulation time: ", time
        write(*,'(a)') "Wrote diagnostics.csv and f_final.dat"

        if (step >= max_steps .and. time < t_end) then
            write(*,'(a)') "WARNING: max_steps reached before t_end."
        end if

        deallocate(x, v, f, fnew, fx, fv, density, rho, electric_field)

    end subroutine run_simulation

end module vlasov_poisson_module


program main
    use vlasov_poisson_module, only: run_simulation
    implicit none

    call run_simulation()
end program main
