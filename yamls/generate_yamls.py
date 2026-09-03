import os

template_filename = "MCMC_1.yaml"
with open(template_filename, "r") as f:
    template = f.read()

for i in range(2, 12):
    new_filename = f"MCMC_{i}.yaml"
    contents = template.replace("datavector_01.modelvector", f"datavector_{i:02d}.modelvector")
    contents = contents.replace("MCMC_1", f"MCMC_{i}")
    with open(new_filename, "w") as f:
        f.write(contents)