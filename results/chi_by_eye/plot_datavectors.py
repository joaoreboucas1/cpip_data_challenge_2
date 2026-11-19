# Plots the 21 challenge cosmic shear data vectors against a fiducial cosmology.
#
# The fiducial model vector is computed with Cosmolike through TEMPLATE.yaml: LCDM-like point (w = -1), NLA IA with
# A1 = 0.6 and eta1 = -0.5, no photo-z or shear calibration offsets, no baryons (Q = 0). Error bars are the square root of the covariance
# diagonal; the chi2 in each title uses the full covariance (cosmic shear part only).
#
# Output: datavectors_fiducial.pdf, one page with all data vectors overlaid, then one page per data vector.
# Each page is an 8x8 grid of source-bin pairs: xi_+ in the upper triangle, xi_- in the lower one.
#
# Usage (cocoa environment, after `source start_cocoa.sh`): python plot_datavectors.py
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from cobaya.yaml import yaml_load_file
from cobaya.model import get_model
import cosmolike_roman_real_interface as ci

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.normpath(os.path.join(SCRIPT_DIR, "../../data"))
LIKELIHOOD = "roman_real.cosmic_shear"
SHEAR_LEN  = 1080
N_BINS     = 8
N_THETA    = 15
N_DVS      = 21
FIDUCIAL = {"As_1e9": 2.1, "ns": 0.96, "H0": 70.0, "omegab": 0.049, "omegam": 0.3, "w": -1.0,
            "roman_A1_1": 0.6, "roman_A1_2": -0.5}  # NLA: A1, eta1
THETA_RANGE = (2.5, 200.0)  # arcmin, bin edges of DC2.dataset
CLIP_SIGMA = 10  # y range of the per-data-vector panels: fiducial +/- CLIP_SIGMA errors (points outside are arrows)

DATA_COLOR  = "#2a78d6"
MODEL_COLOR = "#eb6834"
PAIRS = [(a, b) for a in range(N_BINS) for b in range(a, N_BINS)]

def compute_fiducial():
    os.chdir(os.environ["ROOTDIR"])  # paths in the yaml are relative to the Cocoa root
    info = yaml_load_file(os.path.join(SCRIPT_DIR, "TEMPLATE.yaml"))
    info.pop("output", None)
    # The data vector does not enter the model vector; the local roman_real does not know this option
    info["likelihood"][LIKELIHOOD].pop("data_vector_file", None)
    model = get_model(info)
    point = {p: FIDUCIAL.get(p, 0.0) for p in model.parameterization.sampled_params()}
    model.logposterior(point)
    like = model.likelihood[LIKELIHOOD]
    params = {**point, **model.parameterization.constant_params()}
    dv = np.asarray(like.internal_get_datavector(**params))[:SHEAR_LEN]
    theta = np.asarray(ci.get_binning_real_space())  # arcmin
    return dv, theta

def load_cov():
    table = np.loadtxt(os.path.join(DATA_DIR, "cov_sim.dat"), usecols=(0, 1, 8))
    cov = np.zeros((SHEAR_LEN, SHEAR_LEN))
    i, j = table[:, 0].astype(int), table[:, 1].astype(int)
    cov[i, j] = table[:, 2]
    cov[j, i] = table[:, 2]
    return cov

def segment(vec, pm, a, b):
    """Entries of xi_+ (pm = 0) or xi_- (pm = 1) for source bins (a, b), a <= b"""
    k = PAIRS.index((a, b))
    start = pm*SHEAR_LEN//2 + k*N_THETA
    return vec[start:start + N_THETA]

def make_grid(title):
    fig, axes = plt.subplots(N_BINS, N_BINS, figsize=(16, 16), sharex=True)
    fig.suptitle(title, fontsize=14, y=0.995)
    for r in range(N_BINS):
        for c in range(N_BINS):
            ax = axes[r, c]
            ax.set_xscale("log")
            ax.set_xlim(*THETA_RANGE)  # errorbar would otherwise stretch the log autoscale well below the data
            ax.tick_params(labelsize=6, labelbottom=True)
            ax.set_xlabel(r"$\theta$ [arcmin]", fontsize=7, labelpad=1)
            ax.grid(alpha=0.25, lw=0.5)
            pm = 0 if c >= r else 1
            a, b = (r, c) if c >= r else (c, r)
            ax.text(0.04, 0.9, f"$\\xi_{'+' if pm == 0 else '-'}$ ({a+1},{b+1})", transform=ax.transAxes, fontsize=7)
    fig.text(0.005, 0.5, r"$\theta\,\xi_\pm(\theta)\ [10^{-4}\,\mathrm{arcmin}]$", rotation=90, va="center", fontsize=11)
    return fig, axes

def panels(axes):
    """(ax, pm, a, b) for every panel: xi_+ in the upper triangle (incl. diagonal), xi_- in the lower"""
    for r in range(N_BINS):
        for c in range(N_BINS):
            yield (axes[r, c], 0, r, c) if c >= r else (axes[r, c], 1, c, r)

if __name__ == "__main__":
    model, theta = compute_fiducial()
    cov = load_cov()
    inv_cov = np.linalg.inv(cov)
    sigma = np.sqrt(np.diag(cov))
    dvs = [np.loadtxt(os.path.join(DATA_DIR, f"datavectors/datavector_{i:02d}.modelvector"), usecols=1,
                      max_rows=SHEAR_LEN) for i in range(1, N_DVS + 1)]
    chi2 = [(dv - model) @ inv_cov @ (dv - model) for dv in dvs]
    scale = theta*1e4
    short = {"roman_A1_1": "A1", "roman_A1_2": "eta1"}
    fid = ", ".join(f"{short.get(k, k)} = {v}" for k, v in FIDUCIAL.items())

    out = os.path.join(SCRIPT_DIR, "datavectors_fiducial.pdf")
    with PdfPages(out) as pdf:
        # Page 1: all data vectors overlaid
        fig, axes = make_grid(f"All {N_DVS} data vectors (thin) vs fiducial\n{fid}")
        for ax, pm, a, b in panels(axes):
            m = segment(model, pm, a, b)*scale
            for dv in dvs:
                ax.plot(theta, segment(dv, pm, a, b)*scale, color=DATA_COLOR, lw=0.7, alpha=0.5)
            ax.plot(theta, m, color=MODEL_COLOR, lw=2)
            span = np.abs(m).max()
            ax.set_ylim(-1.5*span, 3*span)  # some data vectors are far off-scale (see their own pages)
        axes[0, 0].plot([], [], color=DATA_COLOR, lw=0.7, label="data vectors")
        axes[0, 0].plot([], [], color=MODEL_COLOR, lw=2, label="fiducial")
        axes[0, 0].legend(fontsize=7, loc="lower left")
        fig.tight_layout(rect=(0.01, 0, 1, 0.99))
        pdf.savefig(fig)
        plt.close(fig)

        # One page per data vector
        for i, dv in enumerate(dvs, start=1):
            fig, axes = make_grid(f"Data vector {i:02d} vs fiducial, $\\chi^2$ = {chi2[i-1]:.4g} for {SHEAR_LEN} points\n"
                                  f"{fid}")
            for ax, pm, a, b in panels(axes):
                m, d, s = segment(model, pm, a, b)*scale, segment(dv, pm, a, b)*scale, segment(sigma, pm, a, b)*scale
                lo, hi = (m - CLIP_SIGMA*s).min(), (m + CLIP_SIGMA*s).max()
                inside = (d >= lo) & (d <= hi)
                ax.plot(theta, m, color=MODEL_COLOR, lw=1.5)
                ax.errorbar(theta[inside], d[inside], yerr=s[inside], fmt="o", ms=2.5, color=DATA_COLOR, lw=0.8)
                ax.plot(theta[d > hi], np.full((d > hi).sum(), hi), "^", ms=5, color=DATA_COLOR)
                ax.plot(theta[d < lo], np.full((d < lo).sum(), lo), "v", ms=5, color=DATA_COLOR)
                ax.set_ylim(lo, hi)
            axes[0, 0].plot([], [], "o", ms=3, color=DATA_COLOR, label=f"DV{i:02d}")
            axes[0, 0].plot([], [], "^", ms=5, color=DATA_COLOR, label=f"off-scale (> {CLIP_SIGMA}$\\sigma$)")
            axes[0, 0].plot([], [], color=MODEL_COLOR, lw=1.5, label="fiducial")
            axes[0, 0].legend(fontsize=6, loc="lower left")
            fig.tight_layout(rect=(0.01, 0, 1, 0.99))
            pdf.savefig(fig)
            plt.close(fig)

    for i, c in enumerate(chi2, start=1):
        print(f"DV{i:02d}: chi2 vs fiducial = {c:.4g}")
    print(f"Saved {out}")
