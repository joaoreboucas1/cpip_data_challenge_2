import os

template_filename = "MCMC_1.yaml"
with open(template_filename, "r") as f:
    template = f.read()

# Challenge runs: one per data vector
for i in range(2, 22):
    new_filename = f"MCMC_{i}.yaml"
    contents = template.replace("datavector_01.modelvector", f"datavector_{i:02d}.modelvector")
    contents = contents.replace("MCMC_1", f"MCMC_{i}")
    with open(new_filename, "w") as f:
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
        f.write(contents)
