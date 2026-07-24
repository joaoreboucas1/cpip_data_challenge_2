# Roman Data Challenge 2
This is my attempt at the Roman DC2.

I am using Cocoa v4.11.2.

## Issues and Feedbacks
- Cosmolike needs as input the full data vector, rather than the cosmic shear only part which was provided. The user needs to fill the other entries manually. While the standard computation for a lens-equals-source setup with 8 tomographic bins and 15 angular bins would yield 2160 elements, the `cocoa_roman_real` project excludes 3 GGL bin pairs by default, and then the expected data vector length is 2115 (= 2160 - 3*15). The solution is to pad the given data vectors with zeroes, see `data/datavectors/fill_datavectors.py`. While Cosmolike crashes when the provided data vector has invalid length, it could be useful to log the expected length and the actual length during the crash.
- Additionally, the covariance matrix is not given in a Cosmolike-supported format. Cosmolike supports covariance files with either 3 or 10 columns. The covariance format itself is fine, however, in its header, it states that: "# Column 9 (non-Gaussian) is omitted: a sample covariance yields only the total, not the Gaussian/non-Gaussian split." A temporary solution would be to append a new column filled with zeros to the given covariance, since Cosmolike adds columns 8 and 9, see `cosmolike/generic_interface.cpp`, function `set_inv_cov`. This is implemented in `data/fill_cov.py`.
- When trying to run a likelihood evaluation with `ones.mask`, Cosmolike raises an error regarding that the matrix is singular:
    ```
    File "/Users/joao/cosmo/cocoa/Cocoa/cobaya/cobaya/likelihoods/roman_real/_cosmolike_prototype_base.py", line 124, in initialize
        ci.init_data_real(self.cov_file, self.mask_file, self.data_vector_file)
    RuntimeError: inv(): matrix is singular
    ```
    However, trying to invert the covariance with Numpy in `data/check_covmat.py` is successful, indicating that the matrix is not singular, but actually we need to fill the 3x2pt part just like with the data vectors. A temporary solution would be to fill the diagonal entries (i, i), with i going from 780 to 2114, with ones. This is also done in `data/fill_cov.py` script.

-----------
# Original README Below
-----------

# Data Challenge 2 - Cosmic shear
In the Roman CPIP Data Challenge 2, we aim to
- Perform blind analyses with Cocoa on realistic medium-tier data vectors
- Validate and compare available methods for dealing with baryonic feedback and intrinsic alignments

For the main task of the challenge, participants will receive 25 real space cosmic shear data vectors, which correspond to different wCDM cosmologies with distinct IA and baryonic feedback scenarios. The data vectors were measured from realistic mocks, some of them being dark-matter-only, and others being baryonified mocks. The participants should correctly infer the cosmology of the full set of data vectors. In order to do this, they are encouraged to test and compare different approaches for baryonic feedback, reporting it also in their findings.


## Survey specifications:
The Data Challenge 2 will prepare for the analysis of the medium-tier HLWAS. The survey footprint area is 2330 square degrees, and the effective number density per redshift bin is 37.78 galaxies per squared arcminutes.


## Data vector specifications:
The data vector has length 1080, the first half being $\xi_+$, and the second half being $\xi_-$. For each correlation function, the data vector includes all auto- and cross-correlations between the 8 redshift bins. For each bin combination, there are 15 values for the correlation functions, corresponding to logarithmically spaced values of $\theta$, with $\theta_{\mathrm{min}}=2.5'$ and $\theta_{\mathrm{max}}=200'$.

## Covariance specifications:
Our covariance matrix was computed from a set of 4560 realizations. The format of the file is:

column 0, 1: The covariance matrix indices

column 2, 3: The mean angle $\theta$ of the element

column 4, 5, 6, 7: The tomographic bin indices

column 8: The total value of the covariance

## Timing:
Our challenge has now gone through a soft start. Before we go through an official start and impose a due date, we would like the participants to start running some chains with the data and providing feedback on whether everything is working.

## Instructions:
- Make sure your cocoa installation is up to date
- Build your yaml files with your preferred analysis choices
- Use the provided data vectors, covariance and n(z) to obtain cosmological constraints.
- Submit your feedback on any issue you find during testing to Rafael ([rchgomes@sas.upenn.edu](mailto:rchgomes@sas.upenn.edu)), Megan ([zhaoyif@sas.upenn.edu](mailto:zhaoyif@sas.upenn.edu)), and Jiachuan ([jiac.xu@northeastern.edu](mailto:jiac.xu@northeastern.edu)).
- After the official start, obtain the constraints for all data vectors and submit your results as a pdf with (1) Constraints for all cosmological parameters and nuisance parameters (you may marginalize over the baryonic feedback parameters); and (2) A short description of your analysis choices and additional tests in case you perform them. Additional results may include, for example, an estimate of the baryon suppression ratio for one of the data vectors, or even attempts to identify which data vectors have baryonic feedback and which ones are dark matter only.
- The results should be submitted to Rafael ([rchgomes@sas.upenn.edu](mailto:rchgomes@sas.upenn.edu)), Megan ([zhaoyif@sas.upenn.edu](mailto:zhaoyif@sas.upenn.edu)), and Jiachuan ([jiac.xu@northeastern.edu](mailto:jiac.xu@northeastern.edu)).


### *Validations of data vector done*
- analyzed the 2pcf for the fiducial mocks✅
- run full chain (to do)
