import numpy as np

covmat = np.loadtxt("cov_sim.dat")

# Adding an extra column
zeros = np.zeros(len(covmat))
new_covmat = np.column_stack((covmat, zeros))

# Adding extra lines for the diagonal entries
new_lines = []

shear_len = 1080
full_len = 2115

for i in range(shear_len, full_len):
    new_lines.append([i, i, 1.0, 1.0, 0, 0, 0, 0, 1.0, 0.0])

new_covmat = np.vstack((new_covmat, np.asarray(new_lines)))

print(new_covmat.shape)

np.savetxt("cov_sim_cosmolike_format.dat", new_covmat, fmt="%d %d %.8e %.8e %d %d %d %d %.8e %.8e")
