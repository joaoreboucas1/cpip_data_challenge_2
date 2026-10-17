# Baryonic PCs from the Antilles simulations, cosmic shear only, noise-weighted.
#
# Same method as compute_pcs_manual.ipynb (sklearn PCA, 4 components, of the data-vector differences antilles-i - DMO,
# with no mean subtracted so that the DMO data vector maps to Q = 0), except that each element is divided by its
# data error sigma_i (diagonal of the covariance) before the PCA, so that the PCs describe the baryonic modes that
# matter relative to the noise rather than the ones with the largest absolute xi_pm. The saved PCs are the whitened
# PCs multiplied back by sigma_i (data units, as Cosmolike expects: dv += sum_i Q_i PC_i). With this convention
# Q_i = v_i . (diff/sigma) is in units of the data noise.
#
# Also reports the unweighted PCs (the notebook's convention) for comparison, the reconstruction chi2 of the Antilles
# simulations with 1-4 PCs, and priors on Q from the distribution of the Antilles projections plus a buffer.
#
# Usage (cocoa environment, needs sklearn), after compute_antilles_dvs_fast.py: python compute_pcs.py
# Outputs: BARYON_PCA_ANTILLES.txt (2115 rows, 4 columns, zero outside cosmic shear), pcs_summary.txt
import os
import sys
import logging
import numpy as np
from sklearn.decomposition import PCA

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.normpath(os.path.join(SCRIPT_DIR, "../../data"))
SUMMARY    = os.path.join(SCRIPT_DIR, "pcs_summary.txt")
SHEAR_LEN  = 1080
FULL_LEN   = 2115
N_PCS      = 4
N_SIMS     = 400
PRIOR_BUFFER = 0.25  # fraction of the Antilles range of Q added on each side of the prior

logger = logging.getLogger("pcs")
logger.setLevel(logging.INFO)
logger.propagate = False
for handler in (logging.StreamHandler(sys.stdout), logging.FileHandler(SUMMARY, mode="w")):
    handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(handler)

def load_dv(path):
    return np.loadtxt(path, usecols=1, max_rows=SHEAR_LEN)

def load_cov(path):
    table = np.loadtxt(path, usecols=(0, 1, 8))
    cov = np.zeros((SHEAR_LEN, SHEAR_LEN))
    i, j = table[:, 0].astype(int), table[:, 1].astype(int)
    cov[i, j] = table[:, 2]
    cov[j, i] = table[:, 2]
    return cov

def pca_components(diffs):
    """PCs (rows) as in the notebook: sklearn PCA fit (centred), components used without the mean"""
    pca = PCA(n_components=N_PCS)
    pca.fit(diffs)
    return pca.components_, pca.explained_variance_ratio_

if __name__ == "__main__":
    dmo = load_dv(os.path.join(SCRIPT_DIR, "NOBARYONS.modelvector"))
    sims = np.array([load_dv(os.path.join(SCRIPT_DIR, f"antilles_dvs/antilles-{i}.modelvector"))
                     for i in range(1, N_SIMS + 1)])
    diffs = sims - dmo
    cov = load_cov(os.path.join(DATA_DIR, "cov_sim.dat"))
    inv_cov = np.linalg.inv(cov)
    sigma = np.sqrt(np.diag(cov))

    chi2_sims = np.einsum("si,ij,sj->s", diffs, inv_cov, diffs)
    logger.info(f"{N_SIMS} Antilles simulations, cosmic shear only ({SHEAR_LEN} entries)")
    logger.info(f"  chi2(antilles - DMO), full covariance: median {np.median(chi2_sims):.1f}, "
                f"min {chi2_sims.min():.1f}, max {chi2_sims.max():.1f}")

    # Noise-weighted PCs (saved) and unweighted PCs (comparison)
    pcs_white, var_w = pca_components(diffs/sigma)
    pcs_w = pcs_white*sigma  # back to data units
    pcs_unw, var_unw = pca_components(diffs)
    logger.info(f"  weighted PCs: explained variance (whitened space) {np.array2string(var_w, precision=4)}")
    logger.info(f"  unweighted PCs: explained variance {np.array2string(var_unw, precision=4)}")
    cosines = np.abs((pcs_w/np.linalg.norm(pcs_w, axis=1)[:, None]) @ pcs_unw.T)
    logger.info("  |cos| between weighted PCs (rows, data units, normalized) and unweighted PCs (columns):")
    for i, row in enumerate(cosines):
        logger.info(f"    PC{i+1}: " + " ".join(f"{c:.3f}" for c in row))

    # Reconstruction of the simulations with the first n weighted PCs
    Q = (diffs/sigma) @ pcs_white.T  # (N_SIMS, N_PCS)
    logger.info("  chi2 of (antilles - DMO - reconstruction), full covariance, using the first n weighted PCs:")
    for n in range(N_PCS + 1):
        res = diffs - Q[:, :n] @ pcs_w[:n]
        chi2 = np.einsum("si,ij,sj->s", res, inv_cov, res)
        logger.info(f"    n = {n}: median {np.median(chi2):9.3f}, 95th percentile {np.percentile(chi2, 95):9.3f}, "
                    f"max {chi2.max():9.3f}")

    # Priors on Q: Antilles range plus a buffer, always including the DMO point Q = 0
    logger.info(f"  Q of the Antilles simulations (units of data noise) and priors (range +/- {PRIOR_BUFFER:.0%}, "
                f"including Q = 0):")
    for i in range(N_PCS):
        lo, hi = Q[:, i].min(), Q[:, i].max()
        span = hi - lo
        prior_lo, prior_hi = min(lo - PRIOR_BUFFER*span, 0.0), max(hi + PRIOR_BUFFER*span, 0.0)
        logger.info(f"    Q{i+1}: mean {Q[:, i].mean():8.3f}, std {Q[:, i].std():7.3f}, "
                    f"range [{lo:8.3f}, {hi:8.3f}] -> prior [{prior_lo:8.3f}, {prior_hi:8.3f}]")

    full = np.zeros((FULL_LEN, N_PCS))
    full[:SHEAR_LEN] = pcs_w.T
    np.savetxt(os.path.join(SCRIPT_DIR, "BARYON_PCA_ANTILLES.txt"), full)
    np.savetxt(os.path.join(SCRIPT_DIR, "antilles_Q.txt"), Q, header="Q1 Q2 Q3 Q4 of antilles-1 ... antilles-400")
    print(f"Saved BARYON_PCA_ANTILLES.txt, antilles_Q.txt and {SUMMARY}")
