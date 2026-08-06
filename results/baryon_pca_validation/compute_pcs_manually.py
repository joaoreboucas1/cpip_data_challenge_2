import os
import sys
import subprocess

rootdir = os.getenv("ROOTDIR")
os.chdir(rootdir)

def run_cobaya(i):
    print("----------")
    print(f"Running antilles-{i}")
    completed_proc = subprocess.run(
        ["cobaya-run", "projects/cpip_data_challenge_2/yamls/DV_ANTILLES.yaml", "-f"],
        capture_output=True
    )
    if completed_proc.returncode != 0:
        print(f"ERROR generating dv for antilles-{i}")
        exit(1)

# "DV_ANTILLES.yaml" is already configured with "antilles-1"
# Run "DV_ANTILLES.yaml" with Cobaya
run_cobaya(1)

for i in range(2, 401):
    # Modify "DV_ANTILLES.yaml" such that it now uses antilles-2 and saves the data vector with another name
    with open("projects/cpip_data_challenge_2/yamls/DV_ANTILLES.yaml", "r") as f:
        contents = f.read()
    assert f"antilles-{i-1}" in contents
    contents = contents.replace(f"antilles-{i-1}", f"antilles-{i}")

    with open("projects/cpip_data_challenge_2/yamls/DV_ANTILLES.yaml", "w") as f:
        f.write(contents)

    run_cobaya(i)
