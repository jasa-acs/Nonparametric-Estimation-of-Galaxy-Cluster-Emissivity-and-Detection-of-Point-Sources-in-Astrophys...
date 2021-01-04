clear all
addpath(strcat(pwd,'/code'))


for nM=0:2
    
    cont=0;
    rsummary=struct('F',[],'f',[],'d',[],'fhat',[],'fdom',[],'fhatfull',[],'fdomfull',[],'dhat',[],'type',[]);
    
    ind=[0 0 0];
    M=2^(nM+7);  %128,256 or 512
    nk=12*2^(2-nM);   %48, 24, 12
	load(strcat('data_and_results/EOsim',num2str(M),'.mat'));
    
    for i=1:4

        if(i==1) 
            type='blockssym'; %cosmoblocks
        elseif(i==2) 
            type='f4sym'; %cosmo2
        elseif(i==3) 
            type='truef'; %cosmo1
        elseif(i==4) 
            type='f4'; %como2 asymmetric
        end

        for ps=0:1
            
            %filename=strcat('data_and_results/full_astrosim_diagM',num2str(M),type,'PS',num2str(ps),'B1W1S1Erealnk',num2str(nk),'forcepos2.mat');
            filename=strcat('data_and_results/full_astrosim_diag_powerlaw_M',num2str(M),type,'PS',num2str(ps),'B1W1S1Erealnk',num2str(nk),'forcepos2.mat');
            load(filename);
			
			%%%%%%STEP 1
			%read all results of QUTlasso from simulations and:
			%	- Obtain result using state of the art (SA) method
			%	- Obtain MSE, and log-MSE for QUTlasso and SA
			%	- Calculate FDR and TPR for point sources
			%	- Save in each file in variable summary
			
				mse=0;msedom=0;logmse=0;logmsedom=0;                                          
                msec=0;msedomc=0;logmsec=0;logmsedomc=0;
                msem=0;msedomm=0;logmsem=0;logmsedomm=0;
                mseb=0;msedomb=0;logmseb=0;logmsedomb=0;
                msefull=0;msedomfull=0;logmsefull=0;logmsedomfull=0;
                
                bias=0;biasdom=0;logbias=0;logbiasdom=0;                                          
                biasc=0;biasdomc=0;logbiasc=0;logbiasdomc=0;
                biasm=0;biasdomm=0;logbiasm=0;logbiasdomm=0;
                biasb=0;biasdomb=0;logbiasb=0;logbiasdomb=0;
                biasfull=0;biasdomfull=0;logbiasfull=0;logbiasdomfull=0;
                
				allfhat=[];allfdom=[];allmse=[];allmsedom=[];alllogmse=[];alllogmsedom=[];allbias=[];allbiasdom=[];alllogbias=[];alllogbiasdom=[];
                allfhatc=[];allfdomc=[];allmsec=[];allmsedomc=[];alllogmsec=[];alllogmsedomc=[];allbiasc=[];allbiasdomc=[];alllogbiasc=[];alllogbiasdomc=[];
                allfhatm=[];allfdomm=[];allmsem=[];allmsedomm=[];alllogmsem=[];alllogmsedomm=[];allbiasm=[];allbiasdomm=[];alllogbiasm=[];alllogbiasdomm=[];
                allfhatb=[];allfdomb=[];allmseb=[];allmsedomb=[];alllogmseb=[];alllogmsedomb=[];allbiasb=[];allbiasdomb=[];alllogbiasb=[];alllogbiasdomb=[];
                allfhatfull=[];allfdomfull=[];allmsefull=[];allmsedomfull=[];alllogmsefull=[];alllogmsedomfull=[];allbiasfull=[];allbiasdomfull=[];alllogbiasfull=[];alllogbiasdomfull=[];
             
				FDR=0;TPR=0;
				blurthresh=0.25;
				bluralpha=1.449;
				blurR0=2.2364;
				
				[K,B,qmf,S]=setUpOperatorsWS(M,M^2,M/2,bluralpha,blurR0,blurthresh,'Haar',0,2,2);
				
				for j=1:nk
				
					%Get simulated Image
					F=result(j).F;
                    
                    xqut=linspace(1,M/2*sqrt(2),M/2);
                    xqutfull=linspace(-M/2*sqrt(2),M/2*sqrt(2),M);
                    
                    %oracle profile
                    f=result(j).f;
                    fi=fliplr(f(1:(M/2)));
					fd=f((M/2+1):end);
                    ffull=f;
                    
					f=(fi+fd)/2;
                    fi=interp1(xqut,fi,1:(M/2),'linear');
                    fd=interp1(xqut,fd,1:(M/2),'linear');
                    f=interp1(xqut,f,1:(M/2),'linear');
                    ffull=interp1(xqutfull,ffull,linspace(-M/2,M/2,M),'linear');
                    
                    intervals=[1 floor(M/6) floor(M/3) floor(2*M/3) floor(5*M/6) M];
                    fc=ffull(intervals(3):intervals(4)); %center
                    fm=[ffull(intervals(2):intervals(3)) ffull(intervals(4):intervals(5))]; %middle
                    fb=[ffull(intervals(1):intervals(2)) ffull(intervals(5):intervals(6))]; %border
                    
					%get single estimate (average of left and right estimates)
					fhat=result(j).sol.fhatpos;
					fhati=fliplr(fhat(1:(M/2))')';
					fhatd=fhat((M/2+1):end);
                    fhatfull=fhat;
					fhat=(fhati+fhatd)/2;
                    
					fhati=interp1(xqut,fhati,1:(M/2),'linear');
                    fhatd=interp1(xqut,fhatd,1:(M/2),'linear');
                    fhat=interp1(xqut,fhat,1:(M/2),'linear');
                    fhatfull=interp1(xqutfull,fhatfull,linspace(-M/2,M/2,M),'linear');
                    
                    fhatc=fhatfull(intervals(3):intervals(4)); %center
                    fhatm=[fhatfull(intervals(2):intervals(3)) fhatfull(intervals(4):intervals(5))]; %middle
                    fhatb=[fhatfull(intervals(1):intervals(2)) fhatfull(intervals(5):intervals(6))]; %border
                    
                    result(j).fhati=fhati; %left
                    result(j).fhatd=fhatd; %right
                    result(j).fhat=fhat; %mean left right
                    result(j).fhatfull=fhatfull; %full -M:M
                    result(j).fhatc=fhatc; %center section
                    result(j).fhatm=fhatm; %middle section
                    result(j).fhatb=fhatb; %border section
                    
             
					%Calculate mask for point sources
					Dps=reshape(result(j).d,M,M);
                    %Dps=(Dps>max(Dps(:)/10));
                    Dps=(Dps~=0);
					Dps=(filter2(B,Dps));
					
					if(max(max(Dps))==0)
						Dps=1-Dps;
					else
						Dps(Dps~=0)=1;
						Dps=1-Dps;
					end
					
					%Obtain state-of-the-art estimate (dom) adding a mask in point sources
					Eps=E.*Dps;
					dominique=astroStateOfTheArt(F,Eps,O*0,[0 360]);
					result(j).dominique=dominique;
					
                    fdom=interp1(dominique.xbin,dominique.fbin,1:(M/2),'linear','extrap')';
                    if sum(abs(fdom)==0)
                        fdom=fdom+1e-100;
                    else
                        fdom(fdom<=0)=min(dominique.fbin(dominique.fbin>0));
                    end
                    
                    fdom=fdom';
					result(j).fdom=fdom;
                    
                    dominique=astroStateOfTheArt(F,Eps,O*0,[0 180]);
					result(j).dominiqued=dominique;                    
					
                    fdomd=interp1(dominique.xbin,dominique.fbin,1:(M/2),'linear','extrap')';
                    if sum(abs(fdomd)==0)
                        fdomd=fdomd+1e-100;
                    else
                        fdomd(fdomd<=0)=min(dominique.fbin(dominique.fbin>0));
                    end

                    fdomd=fdomd';
					result(j).fdomd=fdomd;
                    
                    dominique=astroStateOfTheArt(F,Eps,O*0,[180 360]);
					result(j).dominiquei=dominique;    
                    
                    fdomi=interp1(dominique.xbin,dominique.fbin,1:(M/2),'linear','extrap')';
                    if sum(abs(fdomi)==0)
                        fdomi=fdomi+1e-100;
                    else
                        fdomi(fdomi<=0)=min(dominique.fbin(dominique.fbin>0));
                    end
                    
                    fdomi=fdomi';
					result(j).fdomi=fdomi;
                    
                    fdomfull=[fliplr(fdomi) fdomd];
                    fdomfull=interp1([(-M/2):(-1) 1:(M/2)],fdomfull,linspace(-M/2,M/2,M),'linear','extrap');
                    
                    fdomc=fdomfull(intervals(3):intervals(4)); %center
                    fdomm=[fdomfull(intervals(2):intervals(3)) fdomfull(intervals(4):intervals(5))]; %middle
                    fdomb=[fdomfull(intervals(1):intervals(2)) fdomfull(intervals(5):intervals(6))]; %border
                    
                    result(j).fdomfull=fdomfull;
                    
					%obtain mse and logmse of both methods symmetric/center/middle/border/full
					result(j).mse=sum((fhat-f).^2/(M/2));
                    result(j).msec=sum((fhatc-fc).^2/(length(fc)));
                    result(j).msem=sum((fhatm-fm).^2/(length(fm)));
                    result(j).mseb=sum((fhatb-fb).^2/(length(fb)));
                    result(j).msefull=sum((fhatfull-ffull).^2/M);
                    
					result(j).msedom=sum((fdom-f).^2/(M/2));
                    result(j).msedomc=sum((fdomc-fc).^2/(length(fc)));
                    result(j).msedomm=sum((fdomm-fm).^2/(length(fm)));
                    result(j).msedomb=sum((fdomb-fb).^2/(length(fb)));
                    result(j).msedomfull=sum((fdomfull-ffull).^2/M);
                    
					result(j).logmse=sum((log(fhat)-log(f)).^2/(M/2));
                    result(j).logmsec=sum((log(fhatc)-log(fc)).^2/(length(fc)));
                    result(j).logmsem=sum((log(fhatm)-log(fm)).^2/(length(fm)));
                    result(j).logmseb=sum((log(fhatb)-log(fb)).^2/(length(fb)));
                    result(j).logmsefull=sum((log(fhatfull)-log(ffull)).^2/M);
                    
					result(j).logmsedom=sum((log(fdom)-log(f)).^2/(M/2));
                    result(j).logmsedomc=sum((log(fdomc)-log(fc)).^2/(length(fc)));
                    result(j).logmsedomm=sum((log(fdomm)-log(fm)).^2/(length(fm)));
                    result(j).logmsedomb=sum((log(fdomb)-log(fb)).^2/(length(fb)));
                    result(j).logmsedomfull=sum((log(fdomfull)-log(ffull)).^2/M); 
                    

					mse=mse+result(j).mse;
					msedom=msedom+result(j).msedom;
					logmse=logmse+result(j).logmse;
					logmsedom=logmsedom+result(j).logmsedom;
                    
                    msefull=msefull+result(j).msefull;
					msedomfull=msedom+result(j).msedomfull;
					logmsefull=logmsefull+result(j).logmsefull;
					logmsedomfull=logmsedomfull+result(j).logmsedomfull;
                    
                    msec=msec+result(j).msec;
					msedomc=msedomc+result(j).msedomc;
					logmsec=logmsec+result(j).logmsec;
					logmsedomc=logmsedomc+result(j).logmsedomc;
                    
                    msem=msem+result(j).msem;
					msedomm=msedomm+result(j).msedomm;
					logmsem=logmsem+result(j).logmsem;
					logmsedomm=logmsedomm+result(j).logmsedomm;
                    
                    mseb=mseb+result(j).mseb;
					msedomb=msedomb+result(j).msedomb;
					logmseb=logmseb+result(j).logmseb;
					logmsedomb=logmsedomb+result(j).logmsedomb;
                    
                    %%%%%%%

                    %obtain BIAS and logBIAS of both methods symmetric/center/middle/border/full
					result(j).bias=sum((fhat-f)/(M/2));
                    result(j).biasc=sum((fhatc-fc)/(length(fc)));
                    result(j).biasm=sum((fhatm-fm)/(length(fm)));
                    result(j).biasb=sum((fhatb-fb)/(length(fb)));
                    result(j).biasfull=sum((fhatfull-ffull)/M);
                    
					result(j).biasdom=sum((fdom-f)/(M/2));
                    result(j).biasdomc=sum((fdomc-fc)/(length(fc)));
                    result(j).biasdomm=sum((fdomm-fm)/(length(fm)));
                    result(j).biasdomb=sum((fdomb-fb)/(length(fb)));
                    result(j).biasdomfull=sum((fdomfull-ffull)/M);
                    
					result(j).logbias=sum((log(fhat)-log(f))/(M/2));
                    result(j).logbiasc=sum((log(fhatc)-log(fc))/(length(fc)));
                    result(j).logbiasm=sum((log(fhatm)-log(fm))/(length(fm)));
                    result(j).logbiasb=sum((log(fhatb)-log(fb))/(length(fb)));
                    result(j).logbiasfull=sum((log(fhatfull)-log(ffull))/M);
                    
					result(j).logbiasdom=sum((log(fdom)-log(f))/(M/2));
                    result(j).logbiasdomc=sum((log(fdomc)-log(fc))/(length(fc)));
                    result(j).logbiasdomm=sum((log(fdomm)-log(fm))/(length(fm)));
                    result(j).logbiasdomb=sum((log(fdomb)-log(fb))/(length(fb)));
                    result(j).logbiasdomfull=sum((log(fdomfull)-log(ffull))/M); 
                    

					bias=bias+result(j).bias;
					biasdom=biasdom+result(j).biasdom;
					logbias=logbias+result(j).logbias;
					logbiasdom=logbiasdom+result(j).logbiasdom;
                    
                    biasfull=biasfull+result(j).biasfull;
					biasdomfull=biasdom+result(j).biasdomfull;
					logbiasfull=logbiasfull+result(j).logbiasfull;
					logbiasdomfull=logbiasdomfull+result(j).logbiasdomfull;
                    
                    biasc=biasc+result(j).biasc;
					biasdomc=biasdomc+result(j).biasdomc;
					logbiasc=logbiasc+result(j).logbiasc;
					logbiasdomc=logbiasdomc+result(j).logbiasdomc;
                    
                    biasm=biasm+result(j).biasm;
					biasdomm=biasdomm+result(j).biasdomm;
					logbiasm=logbiasm+result(j).logbiasm;
					logbiasdomm=logbiasdomm+result(j).logbiasdomm;
                    
                    biasb=biasb+result(j).biasb;
					biasdomb=biasdomb+result(j).biasdomb;
					logbiasb=logbiasb+result(j).logbiasb;
					logbiasdomb=logbiasdomb+result(j).logbiasdomb;
                    
                    %%%%%%%
                    
					%obtain fdr and tpr for qutlasso
					ind0=find(result(j).d~=0);
					ind1=find(result(j).sol.dhatpos~=0);
					tp=length(intersect(ind0,ind1));
					fp=length(setdiff(ind1,ind0));
					fn=length(ind0)-tp;

					tpr=tp/(tp+fn)
					fdr=fp/(tp+fp);

					if(isnan(fdr)) fdr=1; end

					%create variables with all information
					result(j).tpr=tpr;
					result(j).fdr=fdr;
					FDR=FDR+fdr;
					TPR=TPR+tpr;
					
                    allfhat=[allfhat;fhat];
                    allfhatc=[allfhatc;fhatc];
                    allfhatm=[allfhatm;fhatm];
                    allfhatb=[allfhatb;fhatb];
                    allfhatfull=[allfhatfull;fhatfull];
                    
                    allfdom=[allfdom;fdom];
                    allfdomc=[allfdomc;fdomc];
                    allfdomm=[allfdomm;fdomm];
                    allfdomb=[allfdomb;fdomb];
                    allfdomfull=[allfdomfull;fdomfull];                    
                    
					allmse=[allmse result(j).mse];
                    allmsec=[allmsec result(j).msec];
                    allmsem=[allmsem result(j).msem];
                    allmseb=[allmseb result(j).mseb];
                    allmsefull=[allmsefull result(j).msefull];
                    
					allmsedom=[allmsedom result(j).msedom];
                    allmsedomc=[allmsedomc result(j).msedomc];
                    allmsedomm=[allmsedomm result(j).msedomm];
                    allmsedomb=[allmsedomb result(j).msedomb];
                    allmsedomfull=[allmsedomfull result(j).msedomfull];
                    
					alllogmse=[alllogmse result(j).logmse];
                    alllogmsec=[alllogmsec result(j).logmsec];
                    alllogmsem=[alllogmsem result(j).logmsem];
                    alllogmseb=[alllogmseb result(j).logmseb];
                    alllogmsefull=[alllogmsefull result(j).logmsefull];
                    
					alllogmsedom=[alllogmsedom result(j).logmsedom];
                    alllogmsedomc=[alllogmsedomc result(j).logmsedomc];
                    alllogmsedomm=[alllogmsedomm result(j).logmsedomm];
                    alllogmsedomb=[alllogmsedomb result(j).logmsedomb];
                    alllogmsedomfull=[alllogmsedomfull result(j).logmsedomfull];
                    
                    allbias=[allbias result(j).bias];
                    allbiasc=[allbiasc result(j).biasc];
                    allbiasm=[allbiasm result(j).biasm];
                    allbiasb=[allbiasb result(j).biasb];
                    allbiasfull=[allbiasfull result(j).biasfull];
                    
					allbiasdom=[allbiasdom result(j).biasdom];
                    allbiasdomc=[allbiasdomc result(j).biasdomc];
                    allbiasdomm=[allbiasdomm result(j).biasdomm];
                    allbiasdomb=[allbiasdomb result(j).biasdomb];
                    allbiasdomfull=[allbiasdomfull result(j).biasdomfull];
                    
					alllogbias=[alllogbias result(j).logbias];
                    alllogbiasc=[alllogbiasc result(j).logbiasc];
                    alllogbiasm=[alllogbiasm result(j).logbiasm];
                    alllogbiasb=[alllogbiasb result(j).logbiasb];
                    alllogbiasfull=[alllogbiasfull result(j).logbiasfull];
                    
					alllogbiasdom=[alllogbiasdom result(j).logbiasdom];
                    alllogbiasdomc=[alllogbiasdomc result(j).logbiasdomc];
                    alllogbiasdomm=[alllogbiasdomm result(j).logbiasdomm];
                    alllogbiasdomb=[alllogbiasdomb result(j).logbiasdomb];
                    alllogbiasdomfull=[alllogbiasdomfull result(j).logbiasdomfull];

				end
				
				%Create summary of all results
				summary.allmse=allmse;
				summary.allmsedom=allmsedom;
				summary.alllogmse=alllogmse;
				summary.alllogmsedom=alllogmsedom;
                
                summary.allmsec=allmsec;
				summary.allmsedomc=allmsedomc;
				summary.alllogmsec=alllogmsec;
				summary.alllogmsedomc=alllogmsedomc;
                
                summary.allmsem=allmsem;
				summary.allmsedomm=allmsedomm;
				summary.alllogmsem=alllogmsem;
				summary.alllogmsedomm=alllogmsedomm;
                
                summary.allmseb=allmseb;
				summary.allmsedomb=allmsedomb;
				summary.alllogmseb=alllogmseb;
				summary.alllogmsedomb=alllogmsedomb;
                
                summary.allmsefull=allmsefull;
				summary.allmsedomfull=allmsedomfull;
				summary.alllogmsefull=alllogmsefull;
				summary.alllogmsedomfull=alllogmsedomfull;

				summary.mse=mse/length(result);
				summary.msedom=msedom/length(result);
				summary.logmse=logmse/length(result);
				summary.logmsedom=logmsedom/length(result);
                
                summary.msec=msec/length(result);
				summary.msedomc=msedomc/length(result);
				summary.logmsec=logmsec/length(result);
				summary.logmsedomc=logmsedomc/length(result);
                
                summary.msem=msem/length(result);
				summary.msedomm=msedomm/length(result);
				summary.logmsem=logmsem/length(result);
				summary.logmsedomm=logmsedomm/length(result);
                
                summary.mseb=mseb/length(result);
				summary.msedomb=msedomb/length(result);
				summary.logmseb=logmseb/length(result);
				summary.logmsedomb=logmsedomb/length(result);
                
                summary.msefull=msefull/length(result);
				summary.msedomfull=msedomfull/length(result);
				summary.logmsefull=logmsefull/length(result);
				summary.logmsedomfull=logmsedomfull/length(result);
                
                summary.allmse=allmse;
				summary.allmsedom=allmsedom;
				summary.alllogmse=alllogmse;
				summary.alllogmsedom=alllogmsedom;
                
                summary.allmsec=allmsec;
				summary.allmsedomc=allmsedomc;
				summary.alllogmsec=alllogmsec;
				summary.alllogmsedomc=alllogmsedomc;
                
                summary.allmsem=allmsem;
				summary.allmsedomm=allmsedomm;
				summary.alllogmsem=alllogmsem;
				summary.alllogmsedomm=alllogmsedomm;
                
                summary.allmseb=allmseb;
				summary.allmsedomb=allmsedomb;
				summary.alllogmseb=alllogmseb;
				summary.alllogmsedomb=alllogmsedomb;
                
                summary.allmsefull=allmsefull;
				summary.allmsedomfull=allmsedomfull;
				summary.alllogmsefull=alllogmsefull;
				summary.alllogmsedomfull=alllogmsedomfull;

				summary.bias=bias/length(result);
				summary.biasdom=biasdom/length(result);
				summary.logbias=logbias/length(result);
				summary.logbiasdom=logbiasdom/length(result);
                
                summary.biasc=biasc/length(result);
				summary.biasdomc=biasdomc/length(result);
				summary.logbiasc=logbiasc/length(result);
				summary.logbiasdomc=logbiasdomc/length(result);
                
                summary.biasm=biasm/length(result);
				summary.biasdomm=biasdomm/length(result);
				summary.logbiasm=logbiasm/length(result);
				summary.logbiasdomm=logbiasdomm/length(result);
                
                summary.biasb=biasb/length(result);
				summary.biasdomb=biasdomb/length(result);
				summary.logbiasb=logbiasb/length(result);
				summary.logbiasdomb=logbiasdomb/length(result);
                
                summary.biasfull=biasfull/length(result);
				summary.biasdomfull=biasdomfull/length(result);
				summary.logbiasfull=logbiasfull/length(result);
				summary.logbiasdomfull=logbiasdomfull/length(result);
                
                summary.allbias=allbias;
				summary.allbiasdom=allbiasdom;
				summary.alllogbias=alllogbias;
				summary.alllogbiasdom=alllogbiasdom;
                
                summary.allbiasc=allbiasc;
				summary.allbiasdomc=allbiasdomc;
				summary.alllogbiasc=alllogbiasc;
				summary.alllogbiasdomc=alllogbiasdomc;
                
                summary.allbiasm=allbiasm;
				summary.allbiasdomm=allbiasdomm;
				summary.alllogbiasm=alllogbiasm;
				summary.alllogbiasdomm=alllogbiasdomm;
                
                summary.allbiasb=allbiasb;
				summary.allbiasdomb=allbiasdomb;
				summary.alllogbiasb=alllogbiasb;
				summary.alllogbiasdomb=alllogbiasdomb;
                
                summary.allbiasfull=allbiasfull;
				summary.allbiasdomfull=allbiasdomfull;
				summary.alllogbiasfull=alllogbiasfull;
				summary.alllogbiasdomfull=alllogbiasdomfull;
                
				summary.FDR=FDR/length(result);
				summary.TPR=TPR/length(result);
                
				summary.allfhat=allfhat;
				summary.allfdom=allfdom;
                summary.allfhatc=allfhatc;
				summary.allfdomc=allfdomc;
                summary.allfhatm=allfhatm;
				summary.allfdomm=allfdomm;
                summary.allfhatb=allfhatb;
				summary.allfdomb=allfdomb;
                summary.allfhatfull=allfhatfull;
				summary.allfdomfull=allfdomfull;
                
                summary.f=f;
                summary.fc=fc;
                summary.fm=fm;
                summary.fb=fb;
                summary.ffull=ffull;

				%save summary
				save(filename,'result','summary');
				
			%%%%%%%%%%%%%
			
			
			%%%%%%STEP 2
			% For each simulation, obtain the best, the median, and the worst scenarios and summarize results
			
				ind(1)=find(summary.alllogmsefull==min(summary.alllogmsefull),1);
				ind(2)=find(abs(median(summary.alllogmsefull)-summary.alllogmsefull)==min(abs(median(summary.alllogmsefull)-summary.alllogmsefull)),1);
				ind(3)=find(summary.alllogmsefull==max(summary.alllogmsefull),1);

				for k=1:3
					rsummary(cont+k).F=result(ind(k)).F;
					rsummary(cont+k).f=result(ind(k)).f;
					rsummary(cont+k).d=result(ind(k)).d;
					rsummary(cont+k).fhat=result(ind(k)).fhat;
					rsummary(cont+k).fdom=result(ind(k)).fdom;
                    rsummary(cont+k).fhatfull=result(ind(k)).fhatfull;
					rsummary(cont+k).fdomfull=result(ind(k)).fdomfull;
					rsummary(cont+k).dhat=result(ind(k)).sol.dhatpos;
					rsummary(cont+k).type=(k-1)/2;
					rsummary(cont+k).ps=(sum(result(ind(k)).d)>0);
				end
				
				cont=cont+3;
			%%%%%%%%%%%%%
            
        end
    end
    save(strcat('data_and_results/rsummary_powerlaw_',num2str(M),'.mat'),'rsummary')
end

