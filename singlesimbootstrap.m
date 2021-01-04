WavePath
addpath(strcat(pwd,'/astro/FISTA'))
addpath(strcat(pwd,'/astro/matlab2016'))

load EandOforsim.mat
options=astro_setparams();
options.showplot=0;
options.parallel=1;

M=512
nboot=48
parpool(6)
[sol, F, f, d] = astrosimPS(1,M,M/2,Esim,0,[1000 10 0.55],'truef',0.002,nboot,0,options);
save(strcat('singlesimM',num2str(M),'nboot',num2str(nboot),'.mat'),'sol','F','f','d')