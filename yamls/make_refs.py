# Writes refs/MCMC_<i>.yaml (ref: loc = best fit) from the minimizer outputs MIN_<i>.minimum.txt
#
# Usage: python make_refs.py <folder with the MIN_<i>.minimum.txt files>
#        then python generate_yamls.py to add the refs to MCMC_<i>.yaml
import os
import sys
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLED = ["As_1e9", "ns", "H0", "omegab", "omegam", "w"] \
    + [f"roman_DZ_S{i}" for i in range(1, 9)] + [f"roman_M{i}" for i in range(1, 9)] \
    + ["roman_A1_1", "roman_A1_2"] + [f"roman_BARYON_Q{i}" for i in range(1, 4)]

if __name__ == "__main__":
    folder = sys.argv[1]
    for i in range(1, 22):
        path = os.path.join(folder, f"MIN_{i}.minimum.txt")
        if not os.path.exists(path):
            print(f"MIN_{i}: no minimum file, skipped")
            continue
        with open(path) as f:
            names = f.readline().lstrip("#").split()
            values = [float(v) for v in f.readline().split()]
        best = dict(zip(names, values))
        with open(os.path.join(SCRIPT_DIR, f"refs/MCMC_{i}.yaml"), "w") as out:
            out.write(f"# Reference point of MCMC_{i}: best fit of MIN_{i} (chi2 = {best.get('chi2', float('nan')):.2f})\n")
            yaml.safe_dump({p: {"ref": {"loc": round(best[p], 6)}} for p in SAMPLED}, out, sort_keys=False)
        print(f"MIN_{i}: chi2 = {best.get('chi2', float('nan')):.2f}, omegam = {best['omegam']:.3f}, w = {best['w']:.3f}")
