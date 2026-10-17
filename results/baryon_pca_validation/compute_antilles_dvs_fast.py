# Antilles data vectors in a single process: a faster alternative to compute_antilles_dvs.py.
#
# All Antilles data vectors share the fiducial cosmology; only the baryonic suppression table changes. Instead of
# calling cobaya-run once per simulation, this script builds the Cobaya model from DV_ANTILLES.yaml once, evaluates
# the fiducial point once (CAMB runs once), and then, for each simulation, switches Cosmolike's contamination table
# and recomputes the data vector with the likelihood's own internal_get_datavector (no CAMB call). The DMO data
# vector is computed in the same process (contamination switched off), so it shares the exact same cosmology and
# binning. The yaml is not modified.
#
# Usage (cocoa environment, after `source start_cocoa.sh`):
#   python compute_antilles_dvs_fast.py            # DMO + antilles-1 ... antilles-400
#   python compute_antilles_dvs_fast.py 1 24       # DMO + a range of simulations (inclusive)
# Outputs: antilles_dvs/antilles-{i}.modelvector and NOBARYONS.modelvector (same format as cobaya-run)
import os
import sys
import time
import numpy as np
from cobaya.yaml import yaml_load_file
from cobaya.model import get_model
import cosmolike_roman_real_interface as ci

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LIKELIHOOD = "roman_real.cosmic_shear"

def save_dv(path, dv):
    np.savetxt(path, np.column_stack([np.arange(len(dv)), dv]), fmt=("%d", "%1.8e"))

if __name__ == "__main__":
    first, last = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1, 400)
    os.chdir(os.environ["ROOTDIR"])  # paths in the yaml are relative to the Cocoa root

    yaml_dir = os.path.join(SCRIPT_DIR, "../../yamls")
    info = yaml_load_file(os.path.join(yaml_dir, "DV_ANTILLES.yaml"))
    point = dict(info.pop("sampler")["evaluate"]["override"])
    info.pop("output", None)
    info["likelihood"][LIKELIHOOD]["print_datavector"] = False
    model = get_model(info)

    like = model.likelihood[LIKELIHOOD]
    missing = [p for p in model.parameterization.sampled_params() if p not in point]
    # The PCA amplitudes are not all in the yaml's override; they do not enter when use_baryon_pca is off
    if not like.use_baryon_pca:
        for p in [p for p in missing if p.startswith("roman_BARYON_Q")]:
            point[p] = 0.0
            missing.remove(p)
    if missing:
        raise SystemExit(f"Parameters without a value in the yaml's evaluate override: {missing}")
    model.logposterior(point)  # computes the cosmology (CAMB) once
    params = {**point, **model.parameterization.constant_params()}
    print("theta bins (arcmin):", np.round(np.asarray(ci.get_binning_real_space()), 3), flush=True)

    # Cosmolike caches its results under a key (cosmology.random) that is only renewed when a cosmological input
    # changes, not when the contamination table does. Setting the distances to a slightly different chi(z) and
    # back renews it twice and leaves the exact distances in place.
    z_1D = like.z_interp_1D
    chi = like.provider.get_comoving_radial_distance(z_1D)*like.provider.get_param("H0")/100  # Mpc/h
    def clear_cosmolike_cache():
        ci.set_distances(z=z_1D, chi=chi*(1 + 1e-9))
        ci.set_distances(z=z_1D, chi=chi)

    # DMO: no contamination
    ci.reset_bary_struct()
    clear_cosmolike_cache()
    save_dv(os.path.join(SCRIPT_DIR, "NOBARYONS.modelvector"), np.asarray(like.internal_get_datavector(**params)))
    print("DMO done", flush=True)

    out_dir = os.path.join(SCRIPT_DIR, "antilles_dvs")
    start = time.time()
    for i in range(first, last + 1):
        clear_cosmolike_cache()
        ci.init_baryons_contamination(sim=f"antilles-{i}", allsims=like.allsims)
        save_dv(os.path.join(out_dir, f"antilles-{i}.modelvector"), np.asarray(like.internal_get_datavector(**params)))
        print(f"antilles-{i} done ({(time.time() - start)/(i - first + 1):.2f} s per data vector)", flush=True)
