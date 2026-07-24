# Roman Data Challenge 2
This is my attempt at the Roman DC2.

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
