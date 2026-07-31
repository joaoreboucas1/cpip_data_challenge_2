r"""
Script to compute scale cuts using DES-Y3 methodology (Krause et al 2021, https://arxiv.org/pdf/2105.13548):
- Cosmic shear cuts are determined by finding a values for \theta_min in each pair of bins (i, j) such that, for every bin pair, \Delta\chi^2_(i,j) in that specific bin it less than \Delta\chi^2_threshold / NUM_SHEAR_2PCFS. The \Delta\chi^2 is the difference between a data vector computed with a fiducial model (e.g. no baryonic feedback) and a data vector computed with a model that has a small-scale modification (e.g. baryonic feedback). \Delta\chi^2_threshold is an user input and NUM_SHEAR_2PCFS is equal to 30 (15 xi_plus and 15 xi_minus) for 5 source bins.
- Clustering and GGL scale cuts are detemined by setting minimum separation scales in Mpc/h, and then converting those to minimum angular separations using the average redshift of each lens bin.
The script prints a mask file for Cocoa analyses, so the user can redirect that to a file of their preference.
Usage:
    python scale_cuts.py dv1 dv2 cov [--chi2_threshold threshold] [--Rmin_gc Rmin_gc] [--Rmin_ggl Rmin_ggl] [--output mask_filename]
Where:
    dv1: Path to first data vector (e.g. linear theory)
    dv2: Path to second data vector (e.g. halofit)
    cov: Path to covariance matrix
    threshold: Chi2 threshold for cosmic shear scale cuts (default: 0.5)
    Rmin_gc: Minimum separation scale for galaxy clustering in Mpc/h (default: 1.0)
    Rmin_ggl: Minimum separation scale for galaxy-galaxy lensing in Mpc/h (default: 1.0)
    mask_filename: Output file for mask (default: mask.txt)
"""

from argparse import ArgumentParser
from astropy.cosmology import FlatLambdaCDM
import os
import numpy as np
import itertools

SHEAR_LEN = 1080
NUM_SRC_BINS = 8
NUM_SHEAR_2PCFS = 72 # 8 bins -> 36 xip and 36 xim

ARCMIN_TO_RAD = 2.90888208665721580e-4
RAD_TO_ARCMIN = 1/ARCMIN_TO_RAD
THETA_MIN_ARCMIN  = 2.5   # Minimum angular scale (in arcminutes)
THETA_MAX_ARCMIN  = 250.  # Maximum angular scale (in arcminutes)
NUM_ANG_BINS = 15         # Number of angular bins

THETA_MIN_RAD = THETA_MIN_ARCMIN * ARCMIN_TO_RAD
THETA_MAX_RAD = THETA_MAX_ARCMIN * ARCMIN_TO_RAD
DLOG_THETA = (np.log(THETA_MAX_RAD) - np.log(THETA_MIN_RAD))/NUM_ANG_BINS
theta = np.zeros(NUM_ANG_BINS)

for i in range(NUM_ANG_BINS):
    THETA_MIN_BIN = np.exp(np.log(THETA_MIN_RAD) + i * DLOG_THETA)
    THETA_MAX_BIN = np.exp(np.log(THETA_MIN_RAD) + (i + 1) * DLOG_THETA)
    theta[i] = (2/3) * (THETA_MAX_BIN**3 - THETA_MIN_BIN**3) / (THETA_MAX_BIN**2 - THETA_MIN_BIN**2)

theta *= RAD_TO_ARCMIN

def make_pairs(max_i):
    # Helper function to map integers to symmetric bin pairs
    # e.g. 0 -> (1, 1); 1 -> (1, 2); 2 -> (1, 3)
    # Input max_i is the number of source bins
    pairs = []
    for i in range(1, max_i+1):
        for j in range(i, max_i+1):
            pairs.append((i, j))
    return pairs

pairs = make_pairs(8)

def invert_cov(cov):
    # Covariance Trick from https://arxiv.org/pdf/2601.00438 Appendix C1
    variances = np.diag(cov)
    Lambda = np.diag(1/np.sqrt(variances))
    matrix_to_invert = Lambda @ cov @ Lambda
    invcov = Lambda @ np.linalg.inv(matrix_to_invert) @ Lambda
    return invcov


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--threshold", required=True, type=float, help="Chi2 threshold for cosmic shear scale cuts")
    parser.add_argument("--output", required=True, type=str, help="Output file for mask")
    parser.add_argument("--verbose", action="store_true", help="Enable tracing of shear scale cuts")
    args = parser.parse_args()

    dv_baseline_filename = "BASELINE.modelvector"
    dv_nobaryons_filename = "STRONGBARYONS.modelvector"
    dv_strongbaryons_filename = "NOBARYONS.modelvector"

    _, dv_baseline      = np.loadtxt(dv_baseline_filename, unpack=True)
    _, dv_nobaryons     = np.loadtxt(dv_nobaryons_filename, unpack=True)
    _, dv_strongbaryons = np.loadtxt(dv_strongbaryons_filename, unpack=True)
    
    print(f"Using data vectors {dv_baseline_filename}, {dv_nobaryons_filename}, and {dv_strongbaryons_filename}. Number of elements = {len(dv_baseline)}")

    # Shear part of data vectors for dchi2 calculations
    dv_baseline      = dv_baseline[:SHEAR_LEN]
    dv_nobaryons     = dv_nobaryons[:SHEAR_LEN]
    dv_strongbaryons = dv_strongbaryons[:SHEAR_LEN]

    # Loading cov in Cosmolike format
    cov_filename = "../cov_sim_cosmolike_format.dat"
    cov = np.zeros((SHEAR_LEN, SHEAR_LEN))

    with open(cov_filename, "r") as f:
        for line in f.read().splitlines():
            if line.startswith("#"): continue
            words = line.split()
            i = int(words[0])
            j = int(words[1])
            if i >= SHEAR_LEN or j >= SHEAR_LEN: continue
            cov[i,j] = float(words[8]) + float(words[9])
            cov[j,i] = cov[i,j]

    print(f"Loaded covariance {cov_filename}")

    print(f"Scale cut settings:")
    print(f"  - chi2_threshold = {args.threshold}")
    print("--------------------")

    # Cosmic shear scale cuts
    delta_strongbaryons = dv_baseline - dv_strongbaryons
    delta_nobaryons = dv_baseline - dv_nobaryons

    invcov = invert_cov(cov)

    def compute_cuts_xi_ij():
        minimum_angles = np.zeros((NUM_SHEAR_2PCFS,))
        full_shear_mask = np.empty(0)
        for i in range(NUM_SHEAR_2PCFS):
            shear_func = "xi_plus" if i < NUM_SHEAR_2PCFS//2 else "xi_minus"
            bin_pair = pairs[i%(NUM_SHEAR_2PCFS//2)]

            mask = np.ones((NUM_ANG_BINS,))
            num_removed_elements = 0
            idx_start = i * NUM_ANG_BINS
            idx_end = (i + 1) * NUM_ANG_BINS

            delta_nobaryons_ij     = delta_nobaryons[idx_start:idx_end]
            delta_strongbaryons_ij = delta_strongbaryons[idx_start:idx_end]

            cov_ij = cov[idx_start:idx_end, idx_start:idx_end].copy()
            invcov_ij = invert_cov(cov_ij)

            chi2_nobaryons_ij     = delta_nobaryons_ij     @ invcov_ij @ delta_nobaryons_ij
            chi2_strongbaryons_ij = delta_strongbaryons_ij @ invcov_ij @ delta_strongbaryons_ij
            chi2_ij = max(chi2_nobaryons_ij, chi2_strongbaryons_ij)

            if args.verbose:
                print("--------------------")
                print(f"Initial chi2 for {shear_func} bin pair {bin_pair} = {chi2_ij}")

            while chi2_ij > args.threshold / NUM_SHEAR_2PCFS:
                mask[num_removed_elements] = 0.0
                num_removed_elements += 1

                kept = np.where(mask > 0)[0]
                
                delta_nobaryons_ij_masked     = delta_nobaryons_ij[kept]
                delta_strongbaryons_ij_masked = delta_strongbaryons_ij[kept]

                cov_ij_masked = cov_ij[np.ix_(kept, kept)]
                invcov_ij_masked = np.linalg.inv(cov_ij_masked)
                
                chi2_nobaryons_ij     = delta_nobaryons_ij_masked     @ invcov_ij_masked @ delta_nobaryons_ij_masked
                chi2_strongbaryons_ij = delta_strongbaryons_ij_masked @ invcov_ij_masked @ delta_strongbaryons_ij_masked
                chi2_ij = max(chi2_nobaryons_ij, chi2_strongbaryons_ij)
                # chi2_ij = chi2_strongbaryons_ij
                
                if args.verbose: 
                    print(f"({bin_pair}) Removed {num_removed_elements} out of {NUM_ANG_BINS} elements. chi2 for {shear_func} bin pair {bin_pair} = {chi2_nobaryons_ij:.4f} (baseline vs no baryons), {chi2_strongbaryons_ij:.4f} (baseline vs strong baryons)")

            full_shear_mask = np.concatenate((full_shear_mask, mask))
            if num_removed_elements < NUM_ANG_BINS: minimum_angles[i] = theta[num_removed_elements] # NOTE: minimum_angles[i] is INCLUDED in the scale cuts
            else: minimum_angles[i] = 999
        return minimum_angles, full_shear_mask

    min_angles, shear_mask = compute_cuts_xi_ij()

    print("Minimum angles for shear 2pcfs:")
    for i, (pair, min_angle) in enumerate(zip(itertools.cycle(pairs), min_angles)):
        shear_func = "xi_plus" if i < NUM_SHEAR_2PCFS//2 else "xi_minus"
        print(f"  - {shear_func} {pair}: {min_angle:.1f}")
    print(f"Original number of elements: {len(shear_mask)}, number of unmasked elements: {len(shear_mask[shear_mask > 0])}")
    
    delta_nobaryons_masked = delta_nobaryons[:NUM_SHEAR_2PCFS*NUM_ANG_BINS]*shear_mask
    delta_strongbaryons_masked = delta_strongbaryons[:NUM_SHEAR_2PCFS*NUM_ANG_BINS]*shear_mask
    invcov_shear = invcov[:NUM_SHEAR_2PCFS*NUM_ANG_BINS, :NUM_SHEAR_2PCFS*NUM_ANG_BINS]
    chi2_nobaryons = delta_nobaryons_masked@invcov_shear@delta_nobaryons_masked
    chi2_strongbaryons = delta_strongbaryons_masked@invcov_shear@delta_strongbaryons_masked

    print(f"Shear chi2: {chi2_nobaryons:.4f} (no baryons), {chi2_strongbaryons:.4f} (strong baryons)")
    
    ROMAN_3X2_DV_LEN = 2115
    full_mask = np.concatenate((shear_mask, np.zeros(ROMAN_3X2_DV_LEN - SHEAR_LEN)))

    np.savetxt(args.output, np.vstack([np.arange(len(full_mask)), full_mask]).T, fmt="%d %d")
    print(f"Saved mask in {args.output}")