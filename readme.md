# Nonparametric Estimation of Galaxy Cluster Emissivity and Detection of Point Sources in Astrophysics With Two Lasso Penalties

# Author Contributions Checklist Form

## Data

### Abstract

The data used in this paper were collected by the XMM-Newton observatory from the European Space Agency. XMM-Newton is a space-bourne telescope that is observing the sky in the X-ray range. The data used in this paper correspond to observations of galaxy clusters through their X-ray emission. 

### Availability

All the datasets used in the paper are publicly available. 

### Description

There are no data rights on the data used in the paper. The data can be retrieved on the XMM-Newton science archive (XSA, http://nxsa.esac.esa.int/nxsa-web/#home). The data have then be processed using the XMM-Newton Science Analysis Software (XMMSAS), which is publicly available (https://www.cosmos.esa.int/web/xmm-newton/download-and-install-sas).

## Code

### Abstract 

The provided code includes functions to reproduce all simulations, real data results, figures and tables in this paper. In particular it includes the main functions to estimate the galaxy cluster’s emissivity and location of point sources from 2D telescope images, and the uncertainty of the estimation.

### Description 

The code is delivered in MATLAB, with MIT license. All codes are available in file ‘ALIAS.rar’. Instructions to download it and all preliminary requirements are available in https://sites.google.com/view/alias-code.

## Instructions for use

### Reproducibility 

To perform all the analysis and generate all figures (and tables) in the paper with the existing simulated results, go to the root folder ‘ALIAS’ in MATLAB and run ‘reproduceAnalysis’. To perform all simulations in the paper, go to the root folder ‘ALIAS’ in MATLAB and run 'reproduceALL'. By default it runs using parallel computations, you can change it by uncommenting options.parallel=0 and replacing 'parfor' by 'for', in all simulation files ('runSimTable1parallel.m', 'runSimFigure3parallel.m', and 'runRealDataParallel.m').

### Replication

To obtain the estimates of a single image use function ‘astrosolveWS.m'. To simulate a telescope image and obtain the estimates use function 'astrosimPS.m'. To get the estimates of a single image and obtain pointwise confidence interval by bootstrap use function 'astrobootWS.m'. 

Real data concerning images from abell2142 (Chandra and XMM telescope)  and from a3667, bullet, and perseus (XMM telescope) is available in file ‘realdata.mat’ as a 1x5 structure containing the telescope image (F), telescope sensitivity (E), and telescope background (O) for each of the 5 real data settings.
