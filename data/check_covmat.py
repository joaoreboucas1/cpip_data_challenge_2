import numpy as np

DV_LEN = 1080
cov = np.zeros((DV_LEN, DV_LEN))

with open("cov_sim.dat", "r") as f:
    lines = f.read().splitlines()

for line in lines:
    if line.startswith("#"): continue
    entries = line.split()
    i = int(entries[0])
    j = int(entries[1])
    cov_ij = float(entries[8])
    cov[i, j] = cov_ij
    cov[j, i] = cov_ij

invcov = np.linalg.inv(cov)