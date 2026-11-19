import os

# Reference points: refs/<RUN>.yaml (only `ref: loc` entries) is appended to the params !defaults list when it exists.
# Validation: fiducial + PC projection of the data vector. Challenge: minimizer best fits (see minimize/).
template_filename = "MCMC_1.yaml"
with open(template_filename, "r") as f:
    template = f.read()
template = template.replace("\n  - refs/MCMC_1\n", "\n")  # template without its own ref file
DEFAULTS = "  - params_source_NLA\n"

def with_ref(contents, run):
    if os.path.exists(f"refs/{run}.yaml"):
        contents = contents.replace(DEFAULTS, DEFAULTS + f"  - refs/{run}\n")
    return contents

# Challenge runs: one per data vector
for i in range(1, 22):
    new_filename = f"MCMC_{i}.yaml"
    contents = template.replace("datavector_01.modelvector", f"datavector_{i:02d}.modelvector")
    contents = contents.replace("MCMC_1", f"MCMC_{i}")
    with open(new_filename, "w") as f:
        f.write(with_ref(contents, f"MCMC_{i}"))

# Minimizers of the challenge runs (best fits -> refs/MCMC_<i>.yaml through make_refs.py)
MINIMIZER = """sampler:
  minimize:
    method: bobyqa
    best_of: 2
"""
sampler_start = template.index("sampler:")
for i in range(1, 22):
    contents = template[:sampler_start] + MINIMIZER
    contents = contents.replace("datavector_01.modelvector", f"datavector_{i:02d}.modelvector")
    contents = contents.replace("MCMC_1", f"MIN_{i}")
    with open(f"MIN_{i}.yaml", "w") as f:
        f.write(contents)

# Baryon PCA validation runs: same model as the challenge runs, fiducial data vectors from DV_*.yaml
validation = {
    1: "./projects/cpip_data_challenge_2/data/scale_cuts/BASELINE.modelvector",       # BAHAMAS T78
    2: "./projects/cpip_data_challenge_2/data/scale_cuts/NOBARYONS.modelvector",      # DMO
    3: "./projects/cpip_data_challenge_2/data/scale_cuts/STRONGBARYONS.modelvector",  # OWLS-AGN T85
}
for i, dv in validation.items():
    new_filename = f"BPCA_VAL_{i}.yaml"
    contents = template.replace("./projects/cpip_data_challenge_2/data/datavectors/datavector_01.modelvector", dv)
    contents = contents.replace("MCMC_1", f"BPCA_VAL_{i}")
    with open(new_filename, "w") as f:
        f.write(with_ref(contents, f"BPCA_VAL_{i}"))
