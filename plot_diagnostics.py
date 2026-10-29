import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("diagnostics.csv", delimiter=",", names=True)

t = data["time"]
field_energy = data["electric_field_energy"]
total_energy = data["total_energy"]
mass_error = data["mass_relative_error"]
energy_error = data["energy_relative_error"]

plt.figure()
plt.semilogy(t, np.maximum(field_energy, 1e-30))
plt.xlabel("Time")
plt.ylabel("Electric-field energy")
plt.title("Vlasov-Poisson: Electric-field energy")
plt.grid(True)
plt.tight_layout()
plt.savefig("field_energy.png", dpi=200)

plt.figure()
plt.plot(t, total_energy)
plt.xlabel("Time")
plt.ylabel("Total energy")
plt.title("Vlasov-Poisson: Total energy")
plt.grid(True)
plt.tight_layout()
plt.savefig("total_energy.png", dpi=200)

plt.figure()
plt.semilogy(t, np.maximum(mass_error, 1e-30), label="Mass error")
plt.semilogy(t, np.maximum(energy_error, 1e-30), label="Energy error")
plt.xlabel("Time")
plt.ylabel("Relative error")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("conservation_errors.png", dpi=200)

print("Saved field_energy.png, total_energy.png, conservation_errors.png")
plt.show()
