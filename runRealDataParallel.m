%Run in real data with different wavelets and bootstrapping
%results are given for the xmm and chandra telescopes

t=tempname();
mkdir(t);
disp(['temporary space for parcluster: ', t]);

% some cleanup
c=parcluster('local');
c.JobStorageLocation=t;
delete(c.Jobs);

WavePath
addpath(strcat(pwd,'/FISTA'))
addpath(strcat(pwd,'/code'))
addpath(strcat(pwd,'/simfunc'))

options=astro_setparams();

simoptions=simoptions-1

	if simoptions<45
		flipwave=0;
	else
		flipwave=1;
		simoptions=simoptions-45;
	end

	ind=mod(simoptions,5)+1;
	simoptions=floor(simoptions/5);
	ww=simoptions+1;

	load('data_and_results/realdata.mat')
	curdata=realdata(ind);
	clear('realdata')

	telescope=curdata.telescope;
	curdata.name

	XMM.R0 = 4.862; %arcsec 
	XMM.alpha = 1.52;
	XMM.pixel = 0.000694444444444445; % degree

	Chandra.R0 = 0.367; % arcsec
	Chandra.alpha = 1.68;
	Chandra.pixel = 0.00027333333333334; %degree

	options.ps=1;
	options.wavelet=ww;  %do with fixed wavelet
	options.remblack=0;
	options.nitr=1000;
	options.showplot=0;
	options.parallel=1;
    %options.parallel=0;
	options.center=0.03;

	if(strcmp(telescope,'chandra'))
		'chandra'
		telparams=Chandra;
	elseif strcmp(telescope,'xmm')
		'xmm'
		telparams=XMM;
	end

	telparams.R0=telparams.R0* 1/3600 * 1/telparams.pixel; %pixels
	options.bluralpha=telparams.alpha;
	options.blurR0=telparams.R0;

	options
	parpool(16)
	if strcmp(curdata.name,'perseus')
		nboot=4
	else
		nboot=6
	end

	if flipwave==0
		%left to right wavelets
		filename=strcat('data_and_results/results_ind',num2str(ind),'boot6_',num2str(ww),'.mat')
		[sol,F]= astrobootWS(curdata.F,256,curdata.E,curdata.O,nboot,options);
		save(filename,'sol','F')
	else
		%right to left wavelets
		filename=strcat('data_and_results/results_ind',num2str(ind),'boot6p_',num2str(ww),'.mat')
		[sol,F]= astrobootWS(curdata.Fp,256,curdata.Ep,curdata.Op,nboot,options);
		save(filename,'sol','F')
	end
	



