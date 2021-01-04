The data and code provided here are supplementary information for the paper 

“Nonparametric estimation of galaxy cluster's emissivity and 
point source detection in astrophysics with two lasso penalties". 

Authors: Jairo Diaz Rodriguez, Dominique Eckert, Hatef Monajemi,
Stephane Paltani and Sylvain Sardy


Abstract of the article: 

Astrophysicists are interested in recovering the 3D gas emissivity
of a galaxy cluster from a 2D image taken by a telescope. A blurring
phenomenon and presence of point sources make this inverse problem even
harder to solve. The current state-of-the-art technique is two step: First identify
the location of potential point sources, then mask these locations and
deproject the data.

We instead model the data as a Poisson generalized linear model (involving
blurring, Abel and wavelets operators) regularized by two lasso penalties
to induce sparse wavelet representation and sparse point sources. The
amount of sparsity is controlled by two quantile universal thresholds. As a
result, our method outperforms the existing one.



PRELIMINARY REQUIREMENTS
- This code is implemented in MATLAB

1. 	Download and extract 'ALIAS.rar', available at http://www.unige.ch/math/folks/sardy/astroRepository.
2. 	Add folders 'ALIAS/FISTA', 'ALIAS/simfunc' and 'ALIAS/code' to MATLAB path.
3. 	Download and Install WAVELAB 850. Available at http://statweb.stanford.edu/~wavelab/
4. 	Compile 'ALIAS/deprojection.C' in your own Operative System and locate it in folder 'ALIAS'. This file 
	allows to do the inversion using the conventional method.


INSTRUCTIONS FOR USAGE
1.	To perform all the analysis and generate all figures (and tables) in the paper with the existing simulated 
	results, go to the root folder ‘ALIAS’ in MATLAB and run ‘reproduceAnalysis’. 
2.	To perform all simulations in the paper, go to the root folder ‘ALIAS’ in MATLAB and run 'reproduceALL'. 
	By default it runs using parallel computations, you can change it by uncommenting options.parallel=0 and 
	replacing 'parfor' by 'for', in all simulation files ('runSimTable1parallel.m', 'runSimFigure3parallel.m', 
	and 'runRealDataParallel.m').
3.  To get the estimates in a single image, use 'astrosolveWS.m'.
4.  To simulate a telescope image and get the estimates, use 'astrosimPS.m'.
5.  To get the estimates of a single image and obtain pointwise confidence interval by bootstrap use 'astrobootWS.m'.

REAL DATA
Real data concerning images from abell2142 (Chandra and XMM telescope)  and from a3667, bullet, and perseus 
(XMM telescope) is available in file ‘realdata.mat’ as a 1x5 structure containing the telescope image (F),
 telescope sensitivity (E), and telescope background (O) for each of the 5 real data settings.


Copyright 2019, Jairo Diaz-Rodriguez (adjairo@uninorte.edu.co)
MIT license
