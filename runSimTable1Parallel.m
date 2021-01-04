WavePath
addpath(strcat(pwd,'/FISTA'))
addpath(strcat(pwd,'/code'))
addpath(strcat(pwd,'/simfunc'))

t=tempname();
mkdir(t);
disp(['temporary space for parcluster: ', t]);

% some cleanup
c=parcluster('local');
c.JobStorageLocation=t;
delete(c.Jobs);

E=1000;
forcepos=2;

simoptions=simoptions-1;
type=mod(simoptions,4);
simoptions=floor(simoptions/4);

nM=mod(simoptions,3);
simoptions=floor(simoptions/3);

PS=mod(simoptions,2);

SW=2;B=1;

M=2^(nM+7);
nk=12*2^(2-nM);


if(M==512) 
   parpool(4);
else
   parpool(12);
end

if(type==0)
    type='blockssym'; %cosmoblocks
    sf=3e-7;
elseif(type==1)
    type='truef'; %cosmo1
    sf=[1000 10 0.55];
elseif(type==2)
    type='f4sym'; %cosmo2
    sf=[10000 1 0.45];
elseif(type==3)
    type='f4'; %cosmo 3 asymmetric
    sf=[9000 1 0.45];
end

if    (SW==0)
    W=-1;   S=1;
elseif(SW==1)
    W=1;    S=0;
elseif(SW==2)
    W=1;    S=1;
end

if(W==0) W=-1; end;




options=astro_setparams();
options.bluralpha=options.bluralpha^B; 
options.wavelet(1)=W;
options.ps=PS;
options.showplot=0;
options.sshrink=S;
options.nitr=1000;
options.forcepos=forcepos;
options.psclean=0;
options.center=0;
%options.parallel=0;
options

load('data_and_results/EandOforsim.mat')
clear Osim
result=struct('sol',[],'F',[],'f',[],'d',[]);

parfor i=1:nk
    i
    [sol, F, f, d] = astrosimPS(i,M,M/2,Esim,0,sf,type,0.002,0,0,options);
    result(i).sol=sol;
    result(i).F=F;
    result(i).f=f;
    result(i).d=d;

end

filename=strcat('astrosim_diag_powerlaw_','M',num2str(M),type,'PS',num2str(PS),'B',num2str(B),'W',num2str(W),'S',num2str(S),'E','real','nk',num2str(nk),'forcepos',num2str(forcepos));
save(strcat('data_and_results/full_',filename,'.mat'),'result');
	






