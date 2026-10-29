import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("f_final.dat", comments="#")

x = data[:, 0]
v = data[:, 1]
f = data[:, 2]

x_unique = np.unique(x)
v_unique = np.unique(v)

F = f.reshape(len(v_unique), len(x_unique))

plt.figure(figsize=(9, 5))
plt.pcolormesh(x_unique, v_unique, F, shading="auto")
plt.xlabel("Position x")
plt.ylabel("Velocity v")
plt.title("Final Vlasov phase-space distribution")
plt.colorbar(label="f(x,v)")
plt.tight_layout()
plt.savefig("phase_space_final.png", dpi=200)
plt.show()
