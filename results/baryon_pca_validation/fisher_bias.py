# Fisher-forecast bias of the baryon PCA validation runs (BPCA_VAL_1/2/3) in the (Omega_m, S8) plane.
#
# Linear (Fisher) approximation around the fiducial point theta_fid of the DV_*.yaml files (wCDM + mead2020, w = -1,
# Q = 0). For each fiducial data vector d (DMO, BAHAMAS T78, OWLS-AGN T85) the posterior mean is predicted as
#     theta_mean = theta_fid + F^-1 D^T C^-1 (d - m(theta_fid)),   F = D^T C^-1 D + F_prior,
# with D the derivatives of the cosmic-shear theory vector m (1080 entries) with respect to all sampled parameters
# (central differences). Gaussian priors enter F_prior exactly; flat priors enter as Gaussians of width
# (max - min)/sqrt(12) so that parameters weakly constrained by shear alone (ns, omegab) stay bounded, as in the MCMC.
# S8 = sigma8 sqrt(Omega_m/0.3) is propagated with its Jacobian (sigma8 from CAMB at each step).
# The 2D bias between two means is b = sqrt(delta^T Sigma_2^-1 delta), delta = difference of the (Omega_m, S8) means
# and Sigma_2 the (Omega_m, S8) block of F^-1. Also reported without the PCs (Q fixed to 0) for comparison.
#
# Usage (cocoa environment, after `source start_cocoa.sh` and the patches in patches/):
#   python fisher_bias.py
# Output: fisher_bias.txt (summary) next to this script
import os
import sys
import numpy as np
from cobaya.yaml import yaml_load_file
from cobaya.model import get_model

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT    = os.path.normpath(os.path.join(SCRIPT_DIR, "../.."))
LIKELIHOOD = "roman_real.cosmic_shear"
SHEAR_LEN  = 1080
OUT        = os.path.join(SCRIPT_DIR, "fisher_bias.txt")

FIDUCIALS = {  # name -> data vector (same files as BPCA_VAL_1/2/3)
    "DMO":      os.path.join(PROJECT, "data/scale_cuts/NOBARYONS.modelvector"),
    "BAHAMAS":  os.path.join(PROJECT, "data/scale_cuts/BASELINE.modelvector"),
    "OWLS-AGN": os.path.join(PROJECT, "data/scale_cuts/STRONGBARYONS.modelvector"),
}

# Finite-difference steps (central differences); the PCA amplitudes enter linearly
STEPS = {"As_1e9": 0.02, "ns": 0.005, "H0": 0.5, "omegab": 0.001, "omegam": 0.003, "w": 0.01,
         **{f"roman_DZ_S{i}": 0.0005 for i in range(1, 9)}, **{f"roman_M{i}": 0.001 for i in range(1, 9)},
         "roman_A1_1": 0.05, "roman_A1_2": 0.05, **{f"roman_BARYON_Q{i}": 1.0 for i in range(1, 5)}}

lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(msg)

def load_cov(path):
    table = np.loadtxt(path, usecols=(0, 1, 8))
    cov = np.zeros((SHEAR_LEN, SHEAR_LEN))
    i, j = table[:, 0].astype(int), table[:, 1].astype(int)
    cov[i, j] = table[:, 2]
    cov[j, i] = table[:, 2]
    return cov

def bias_2d(delta, sigma2):
    return float(np.sqrt(delta @ np.linalg.solve(sigma2, delta)))

if __name__ == "__main__":
    os.chdir(os.environ["ROOTDIR"])  # paths in the yamls are relative to the Cocoa root
    yaml_dir = os.path.join(PROJECT, "yamls")

    # Fiducial point: the evaluate override of the DV_*.yaml files, with Q = 0
    fid = dict(yaml_load_file(os.path.join(yaml_dir, "DV_NOBARYONS.yaml"))["sampler"]["evaluate"]["override"])

    info = yaml_load_file(os.path.join(yaml_dir, "BPCA_VAL_2.yaml"))  # DMO data; the derivatives do not use the data
    info.pop("sampler"); info.pop("output", None)
    model = get_model(info)
    like = model.likelihood[LIKELIHOOD]
    if not like.use_baryon_pca:
        raise SystemExit("use_baryon_pca is off in BPCA_VAL_2.yaml")
    sampled = list(model.parameterization.sampled_params())
    for p in sampled:
        fid.setdefault(p, 0.0)  # roman_BARYON_Q*
    missing = [p for p in sampled if p not in fid]
    if missing:
        raise SystemExit(f"No fiducial value for {missing}")
    constants = model.parameterization.constant_params()
    derived_names = list(model.parameterization.derived_params())

    def evaluate(point):
        """Shear theory vector (1080) and S8 at a point"""
        res = model.logposterior(point, cached=False)
        derived = dict(zip(derived_names, res.derived))
        dv = np.asarray(like.get_datavector(**{**point, **constants}))[:SHEAR_LEN]
        return dv, derived["sigma8"]*np.sqrt(point["omegam"]/0.3)

    cov = load_cov(os.path.join(PROJECT, "data/cov_sim.dat"))
    inv_cov = np.linalg.inv(cov)

    m_fid, s8_fid = evaluate(fid)
    data = {k: np.loadtxt(v, usecols=1, max_rows=SHEAR_LEN) for k, v in FIDUCIALS.items()}
    log(f"Fiducial: {', '.join(f'{p} = {fid[p]}' for p in sampled)}; S8 = {s8_fid:.4f}")
    log("Check, chi2(fiducial data vector - model at fiducial, Q = 0): "
        + ", ".join(f"{k} {(d - m_fid) @ inv_cov @ (d - m_fid):.3g}" for k, d in data.items()))
    # Check that the yaml's data_vector_file override is used (BPCA_VAL_2 = DMO -> chi2 ~ 0)
    full_dv = np.asarray(like.get_datavector(**{**fid, **constants}))
    log(f"Check, chi2 computed by the likelihood for BPCA_VAL_2 (DMO data) at the fiducial: "
        f"{-2*like.compute_logp(full_dv):.3g}")

    # Prior Fisher matrix
    prior_var = []
    params_info = model.info()["params"]
    for p in sampled:
        prior = params_info[p]["prior"]
        if "dist" in prior and prior["dist"] == "norm":
            prior_var.append(prior["scale"]**2)
        else:
            prior_var.append((prior["max"] - prior["min"])**2/12)
    F_prior = np.diag(1/np.array(prior_var))

    def derivatives(scale):
        D = np.zeros((SHEAR_LEN, len(sampled)))
        dS8 = np.zeros(len(sampled))
        for k, p in enumerate(sampled):
            h = STEPS[p]*scale
            plus, minus = evaluate({**fid, p: fid[p] + h}), evaluate({**fid, p: fid[p] - h})
            D[:, k] = (plus[0] - minus[0])/(2*h)
            dS8[k] = (plus[1] - minus[1])/(2*h)
        return D, dS8

    def biases(D, dS8, keep):
        """Mean shifts in (Omega_m, S8), their 2D covariance and the biases, using parameters `keep`"""
        idx = [sampled.index(p) for p in keep]
        Dk = D[:, idx]
        F = Dk.T @ inv_cov @ Dk + F_prior[np.ix_(idx, idx)]
        Finv = np.linalg.inv(F)
        J = np.zeros((2, len(idx)))
        J[0, keep.index("omegam")] = 1.0
        J[1] = dS8[idx]
        sigma2 = J @ Finv @ J.T
        means = {k: J @ (Finv @ (Dk.T @ inv_cov @ (d - m_fid))) for k, d in data.items()}
        return means, sigma2

    no_pca = [p for p in sampled if not p.startswith("roman_BARYON_Q")]
    results = {}
    for scale in (1.0, 0.5):
        log(f"\nComputing derivatives (step scale {scale}) ...")
        D, dS8 = derivatives(scale)
        results[scale] = {"PCA": biases(D, dS8, sampled), "no PCA": biases(D, dS8, no_pca)}

    pairs = [("BAHAMAS", "DMO"), ("OWLS-AGN", "DMO"), ("BAHAMAS", "OWLS-AGN")]
    for label in ("PCA", "no PCA"):
        means, sigma2 = results[1.0][label]
        means_half, sigma2_half = results[0.5][label]
        sig = np.sqrt(np.diag(sigma2))
        log(f"\n=== {label} ({'Q1..Q3 marginalized' if label == 'PCA' else 'Q fixed to 0'}) ===")
        log(f"  Fisher errors: sigma(Omega_m) = {sig[0]:.4f}, sigma(S8) = {sig[1]:.4f}, "
            f"corr = {sigma2[0, 1]/sig[0]/sig[1]:.3f}")
        log("  Predicted mean shifts with respect to the truth (in units of sigma), and 2D bias:")
        for k in data:
            dm = means[k]
            log(f"    {k:9s} dOmega_m = {dm[0]/sig[0]:+7.3f} sigma, dS8 = {dm[1]/sig[1]:+7.3f} sigma, "
                f"2D bias = {bias_2d(dm, sigma2):7.3f}   (half step: {bias_2d(means_half[k], sigma2_half):7.3f})")
        log("  2D bias between the means of the validation runs (difference of means):")
        for a, b in pairs:
            log(f"    {a:8s} - {b:8s}: {bias_2d(means[a] - means[b], sigma2):7.3f}   "
                f"(half step: {bias_2d(means_half[a] - means_half[b], sigma2_half):7.3f})")

    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nSaved {OUT}")
