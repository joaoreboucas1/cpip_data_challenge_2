"""
    Issue: Cosmolike does not like a "shear-only" data vector, it requires passing the whole 3x2pt data vector
    Idea: Fill GGL and GC entries with zeros
    Counting DV entries:
        8 source bins + 8 lens bins, 15 angular bins
        Cosmic Shear: 2*T_8 = 2*8*(8+1)/2 = 72 correlation functions, times 15 entries each = 1080 entries (correct from DV length)
        GGL: 8*8 = 64 correlation functions, times 15 entries each = 960 entries
        GC: 8 correlation functions, times 15 entries each = 120 entries
    Caveat: for some unknown reason, the default settings in Cocoa exclude three GGL correlation functions (pairs L6-S0, L7-S0, L7-S1, see cocoa_roman_real/likelihood/cosmic_shear.yaml). This makes the expected data vector length to be 2160 - 3*15 = 2115
    Solution: add another 1080 entries with zeroes
"""

import os
from pathlib import Path

lines_to_add = []
shear_len = 1080
full_len = 2115
for i in range(shear_len, full_len):
    lines_to_add.append(f"{i} 0.0\n")

datavectors_folder = Path(__file__).parent
for filename in os.listdir(datavectors_folder):
    if not filename.endswith(".modelvector"): continue

    with open(datavectors_folder/filename, "r") as f:
        lines = f.readlines()

    lines += lines_to_add
    
    with open(datavectors_folder/filename, "w") as f:
        f.writelines(lines)