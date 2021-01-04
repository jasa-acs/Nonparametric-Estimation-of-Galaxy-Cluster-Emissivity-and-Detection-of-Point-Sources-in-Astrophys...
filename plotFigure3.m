subplot(2,4,1)
M=256
%axis([XMIN XMAX YMIN YMAX])
%x=linspace(-M/2,M/2,M);
x=linspace(0,1,M/2);
load('data_and_results/rsummary_diag_powerlaw_256blockssymps1type1nboot48.mat')
%plotwithci(summary.fmedianfull,summary.fminfull,summary.fmaxfull,summary.fhatfull,summary.ffull,curresult.fdomfull)
plotwithci(x,summary.fmedian,summary.fmin,summary.fmax,summary.fhat,summary.f,curresult.fdom,summary.fdommedian,summary.fdommin,summary.fdommax)
axis1=axis;
title('cosmoBlocks')
ylabel('N=256')

subplot(2,4,2)
load('data_and_results/rsummary_diag_powerlaw_256truefps1type1nboot48.mat')
%plotwithci(summary.fmedianfull,summary.fminfull,summary.fmaxfull,summary.fhatfull,summary.ffull,curresult.fdomfull)
plotwithci(x,summary.fmedian,summary.fmin,summary.fmax,summary.fhat,summary.f,curresult.fdom,summary.fdommedian,summary.fdommin,summary.fdommax)
axis2=axis;
title('cosmo1')

subplot(2,4,3)
load('data_and_results/rsummary_diag_powerlaw_256f4symps1type1nboot48.mat')
%plotwithci(summary.fmedianfull,summary.fminfull,summary.fmaxfull,summary.fhatfull,summary.ffull,curresult.fdomfull)
plotwithci(x,summary.fmedian,summary.fmin,summary.fmax,summary.fhat,summary.f,curresult.fdom,summary.fdommedian,summary.fdommin,summary.fdommax)
axis3=axis;
title('cosmo2')

subplot(2,4,4)
x=linspace(-1,1,M);
load('data_and_results/rsummary_diag_powerlaw_256f4ps1type1nboot48.mat')
plotwithci(x,summary.fmedianfull,summary.fminfull,summary.fmaxfull,summary.fhatfull,summary.ffull,curresult.fdomfull,summary.fdommedianfull,summary.fdomminfull,summary.fdommaxfull)
axis4=axis;
title('cosmoAsym')

M=512
subplot(2,4,5)
x=linspace(0,1,M/2);
load('data_and_results/rsummary_diag_powerlaw_512blockssymps1type1nboot24.mat')
%plotwithci(summary.fmedianfull,summary.fminfull,summary.fmaxfull,summary.fhatfull,summary.ffull,curresult.fdomfull)
plotwithci(x,summary.fmedian,summary.fmin,summary.fmax,summary.fhat,summary.f,curresult.fdom,summary.fdommedian,summary.fdommin,summary.fdommax)
axis(axis1)
ylabel('N=512')
xlabel('r')

subplot(2,4,6)
load('data_and_results/rsummary_diag_powerlaw_512truefps1type1nboot24.mat')
%plotwithci(summary.fmedianfull,summary.fminfull,summary.fmaxfull,summary.fhatfull,summary.ffull,curresult.fdomfull)
plotwithci(x,summary.fmedian,summary.fmin,summary.fmax,summary.fhat,summary.f,curresult.fdom,summary.fdommedian,summary.fdommin,summary.fdommax)
axis(axis2)
xlabel('r')

subplot(2,4,7)
load('data_and_results/rsummary_diag_powerlaw_512f4symps1type1nboot24.mat')
%plotwithci(summary.fmedianfull,summary.fminfull,summary.fmaxfull,summary.fhatfull,summary.ffull,curresult.fdomfull)
plotwithci(x,summary.fmedian,summary.fmin,summary.fmax,summary.fhat,summary.f,curresult.fdom,summary.fdommedian,summary.fdommin,summary.fdommax)
axis(axis3)
xlabel('r')

subplot(2,4,8)
x=linspace(-1,1,M);
load('data_and_results/rsummary_diag_powerlaw_512f4ps1type1nboot24.mat')
plotwithci(x,summary.fmedianfull,summary.fminfull,summary.fmaxfull,summary.fhatfull,summary.ffull,curresult.fdomfull,summary.fdommedianfull,summary.fdomminfull,summary.fdommaxfull)
axis(axis4)
xlabel('r')