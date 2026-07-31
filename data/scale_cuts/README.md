# Scale Cuts

In this folder I perform the scale cuts necessary to mitigate baryonic effects.

## Procedure

- Generate three data vectors: (i) fiducial baryonic feedback (BAHAMAS T78); (ii) no baryonic feedback; (iii) strong baryonic feedback (OWLS-AGN T87)
- For each pair of bins, set the minimum angle $\theta^{ij}_\mathrm{min}$ such that the maximum delta chi2 between baseline and no baryons or baseline and strong baryons falls below a certain threshold
- Validate the scale cuts by running LCDM chains and asserting that the biases are less than a certain threshold.

## Variations
- DES-Y3 only used delta chi2 between baseline and strong baryons
- Junzhou et al used a fixed threshold per bin pair