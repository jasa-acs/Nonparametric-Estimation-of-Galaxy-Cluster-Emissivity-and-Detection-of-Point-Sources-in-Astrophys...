


for nM=2
	M=2^(nM+7);
    
    load(strcat('EOsim',num2str(M),'.mat'))
    blurthresh=0.25;
    bluralpha=1.449;
    blurR0=2.2364;
    [K,B,qmf,S]=setUpOperatorsWS(M,M^2,M/2,bluralpha,blurR0,blurthresh,'Haar',0,2,2);
    
    if(M==256) 
        nboot=48;
    elseif(M==512)
        nboot=24;
    end
	for ntype=1:4
        if(ntype==1) 
			type='blockssym'; %cosmoblocks
		elseif(ntype==2) 
			type='f4sym'; %cosmo2
		elseif(ntype==3) 
			type='truef'; %cosmo1
        elseif(ntype==4) 
			type='f4'; %cosmo1
        end
        for ps=0:1
            for curtype=1
                filename=strcat('data_and_results/rsummary_diag_powerlaw_',num2str(M),type,'ps',num2str(ps),'type',num2str(curtype),'nboot',num2str(nboot),'.mat');
                load(filename)
                
                %Calculate mask for point sources
                Dps=reshape(curresult.d,M,M);
                Dps=(Dps~=0);
                Dps=(filter2(B,Dps));
                if(max(max(Dps))==0)
                    Dps=1-Dps;
                else
                    Dps(Dps~=0)=1;
                    Dps=1-Dps;
                end
                Eps=E.*Dps;
                
                %bootstrap for state-of-the-art
                [bigFdom,bigFdomfull]=image_bootstrapSA(curresult.F,E,O*0,M,0,4,nboot);
                fdomfull=curresult.fdomfull;
                fdom=curresult.fdom;
                
                M=length(sol.fhatpos);
                intervalsfull=[1 floor(M/6) floor(M/3) floor(2*M/3) floor(5*M/6) M];
                intervals= [1 floor(M/6) floor(M/3) floor(M/2)];
                
                xqutfull=linspace(-M/2*sqrt(2),M/2*sqrt(2),M);
                xfull=linspace(-M/2,M/2,M);
                
                xqut=xqutfull((M/2+1):M);
                x=xfull((M/2+1):M);
                %xqut=linspace(1,M/2*sqrt(2),M/2);
                
                f=curresult.f;
                fi=fliplr(f(1:(M/2)));
				fd=f((M/2+1):end);
                ffull=f;
                
                ffull=interp1(xqutfull,ffull,xfull,'linear');
                fi=fliplr(ffull(1:(M/2)));
				fd=ffull((M/2+1):end);
                f=(fi+fd)/2;
                
                fhatfull=curresult.fhatfull;
                fhati=fliplr(fhatfull(1:(M/2)));
				fhatd=fhatfull((M/2+1):end);
                fhat=(fhati+fhatd)/2;
                %fhat=interp1(xqut,fhat,1:(M/2),'linear');
                %fhatfull=interp1(xqutfull,fhatfull,linspace(-M/2,M/2,M),'linear');
                
                
                %asymmetric
                ffullc=ffull(intervalsfull(3):intervalsfull(4)); %center
                ffullm=[ffull(intervalsfull(2):intervalsfull(3)) ffull(intervalsfull(4):intervalsfull(5))]; %middle
                ffullb=[ffull(intervalsfull(1):intervalsfull(2)) ffull(intervalsfull(5):intervalsfull(6))]; %border
                
                fhatfullc=fhatfull(intervalsfull(3):intervalsfull(4)); %center
                fhatfullm=[fhatfull(intervalsfull(2):intervalsfull(3)) fhatfull(intervalsfull(4):intervalsfull(5))]; %middle
                fhatfullb=[fhatfull(intervalsfull(1):intervalsfull(2)) fhatfull(intervalsfull(5):intervalsfull(6))]; %border
                
                fdomfullc=fdomfull(intervalsfull(3):intervalsfull(4)); %center
                fdomfullm=[fdomfull(intervalsfull(2):intervalsfull(3)) fdomfull(intervalsfull(4):intervalsfull(5))]; %middle
                fdomfullb=[fdomfull(intervalsfull(1):intervalsfull(2)) fdomfull(intervalsfull(5):intervalsfull(6))]; %border
                
                
                
                %symmetric
                fb=f(intervals(3):intervals(4)); %border
                fm=f(intervals(2):intervals(3)); %middle
                fc=f(intervals(1):intervals(2)); %center
                
                fhatb=fhat(intervals(3):intervals(4)); %border
                fhatm=fhat(intervals(2):intervals(3)); %middle
                fhatc=fhat(intervals(1):intervals(2)); %center
                
                fdomb=fdom(intervals(3):intervals(4)); %border
                fdomm=fdom(intervals(2):intervals(3)); %middle
                fdomc=fdom(intervals(1):intervals(2)); %center
                
                
                
                
                bigFfull=sol.Fhatpos;
                bigF=zeros(size(bigFfull,1)*2,M/2);
                bigFdom=bigF;
                
                for jj=1:size(bigFfull,1)
                    curffull=bigFfull(jj,:);
                    curfdomfull=bigFdomfull(jj,:);

                    curffull=interp1(xqutfull,curffull,xfull,'linear');
                    curfi=fliplr(curffull(1:(M/2)));
                    curfd=curffull((M/2+1):end);
                    bigFfull(jj,:)=curffull;
                    bigF((2*jj-1):(2*jj),:)=[curfi;curfd];
                    
                    curfi=fliplr(curfdomfull(1:(M/2)));
                    curfd=curfdomfull((M/2+1):end);
                    bigFdom((2*jj-1):(2*jj),:)=[curfi;curfd];
                end
                
                %%%%QUT%%%%%
                %Asymmetric
                fmedianfull=median(bigFfull);
                fmeanfull=mean(bigFfull);
                fmaxfull=quantile(bigFfull,0.975);
                fminfull=quantile(bigFfull,0.025);
                
                fmedianfullc=fmedianfull(intervalsfull(3):intervalsfull(4)); %center
                fmedianfullm=[fmedianfull(intervalsfull(2):intervalsfull(3)) fmedianfull(intervalsfull(4):intervalsfull(5))]; %middle
                fmedianfullb=[fmedianfull(intervalsfull(1):intervalsfull(2)) fmedianfull(intervalsfull(5):intervalsfull(6))]; %border
                
                fmeanfullc=fmeanfull(intervalsfull(3):intervalsfull(4)); %center
                fmeanfullm=[fmeanfull(intervalsfull(2):intervalsfull(3)) fmeanfull(intervalsfull(4):intervalsfull(5))]; %middle
                fmeanfullb=[fmeanfull(intervalsfull(1):intervalsfull(2)) fmeanfull(intervalsfull(5):intervalsfull(6))]; %border
                
                fmaxfullc=fmaxfull(intervalsfull(3):intervalsfull(4)); %center
                fmaxfullm=[fmaxfull(intervalsfull(2):intervalsfull(3)) fmaxfull(intervalsfull(4):intervalsfull(5))]; %middle
                fmaxfullb=[fmaxfull(intervalsfull(1):intervalsfull(2)) fmaxfull(intervalsfull(5):intervalsfull(6))]; %border
                
                fminfullc=fminfull(intervalsfull(3):intervalsfull(4)); %center
                fminfullm=[fminfull(intervalsfull(2):intervalsfull(3)) fminfull(intervalsfull(4):intervalsfull(5))]; %middle
                fminfullb=[fminfull(intervalsfull(1):intervalsfull(2)) fminfull(intervalsfull(5):intervalsfull(6))]; %border
                
                coveragefull=mean(ffull>=fminfull&ffull<=fmaxfull)*100;
                coveragefullc=mean(ffullc>=fminfullc&ffullc<fmaxfullc)*100;
                coveragefullm=mean(ffullm>=fminfullm&ffullm<fmaxfullm)*100;
                coveragefullb=mean(ffullb>=fminfullb&ffullb<fmaxfullb)*100;
                
                %symmetric
                fmedian=median(bigF);
                fmean=mean(bigF);
                fmax=quantile(bigF,1);
                fmin=quantile(bigF,0);
                
                fmedianb=fmedian(intervals(3):intervals(4)); %border
                fmedianm=fmedian(intervals(2):intervals(3)); %middle
                fmedianc=fmedian(intervals(1):intervals(2)); %center
                
                fmeanb=fmean(intervals(3):intervals(4)); %border
                fmeanm=fmean(intervals(2):intervals(3)); %middle
                fmeanc=fmean(intervals(1):intervals(2)); %center
                
                fmaxb=fmax(intervals(3):intervals(4)); %border
                fmaxm=fmax(intervals(2):intervals(3)); %middle
                fmaxc=fmax(intervals(1):intervals(2)); %center
                
                fminb=fmin(intervals(3):intervals(4)); %border
                fminm=fmin(intervals(2):intervals(3)); %middle
                fminc=fmin(intervals(1):intervals(2)); %center
               
                coverage=mean(f>=fmin&f<=fmax)*100;
                coveragec=mean(fc>=fminc&fc<fmaxc)*100;
                coveragem=mean(fm>=fminm&fm<fmaxm)*100;
                coverageb=mean(fb>=fminb&fb<fmaxb)*100;
         
                %%%%SA%%%%%
                %Asymmetric
                fdommedianfull=median(bigFdomfull);
                fdommeanfull=mean(bigFdomfull);
                fdommaxfull=quantile(bigFdomfull,0.975);
                fdomminfull=quantile(bigFdomfull,0.025);
                
                fdommedianfullc=fdommedianfull(intervalsfull(3):intervalsfull(4)); %center
                fdommedianfullm=[fdommedianfull(intervalsfull(2):intervalsfull(3)) fdommedianfull(intervalsfull(4):intervalsfull(5))]; %middle
                fdommedianfullb=[fdommedianfull(intervalsfull(1):intervalsfull(2)) fdommedianfull(intervalsfull(5):intervalsfull(6))]; %border
                
                fdommeanfullc=fdommeanfull(intervalsfull(3):intervalsfull(4)); %center
                fdommeanfullm=[fdommeanfull(intervalsfull(2):intervalsfull(3)) fdommeanfull(intervalsfull(4):intervalsfull(5))]; %middle
                fdommeanfullb=[fdommeanfull(intervalsfull(1):intervalsfull(2)) fdommeanfull(intervalsfull(5):intervalsfull(6))]; %border
                
                fdommaxfullc=fdommaxfull(intervalsfull(3):intervalsfull(4)); %center
                fdommaxfullm=[fdommaxfull(intervalsfull(2):intervalsfull(3)) fdommaxfull(intervalsfull(4):intervalsfull(5))]; %middle
                fdommaxfullb=[fdommaxfull(intervalsfull(1):intervalsfull(2)) fdommaxfull(intervalsfull(5):intervalsfull(6))]; %border
                
                fdomminfullc=fdomminfull(intervalsfull(3):intervalsfull(4)); %center
                fdomminfullm=[fdomminfull(intervalsfull(2):intervalsfull(3)) fdomminfull(intervalsfull(4):intervalsfull(5))]; %middle
                fdomminfullb=[fdomminfull(intervalsfull(1):intervalsfull(2)) fdomminfull(intervalsfull(5):intervalsfull(6))]; %border
                
                coveragedomfull=mean(ffull>=fdomminfull&ffull<=fdommaxfull)*100;
                coveragedomfullc=mean(ffullc>=fdomminfullc&ffullc<fdommaxfullc)*100;
                coveragedomfullm=mean(ffullm>=fdomminfullm&ffullm<fdommaxfullm)*100;
                coveragedomfullb=mean(ffullb>=fdomminfullb&ffullb<fdommaxfullb)*100;
                
                %symmetric
                fdommedian=median(bigFdom);
                fdommean=mean(bigFdom);
                fdommax=quantile(bigFdom,1);
                fdommin=quantile(bigFdom,0);
                
                fdommedianb=fdommedian(intervals(3):intervals(4)); %border
                fdommedianm=fdommedian(intervals(2):intervals(3)); %middle
                fdommedianc=fdommedian(intervals(1):intervals(2)); %center
                
                fdommeanb=fdommean(intervals(3):intervals(4)); %border
                fdommeanm=fdommean(intervals(2):intervals(3)); %middle
                fdommeanc=fdommean(intervals(1):intervals(2)); %center
                
                fdommaxb=fdommax(intervals(3):intervals(4)); %border
                fdommaxm=fdommax(intervals(2):intervals(3)); %middle
                fdommaxc=fdommax(intervals(1):intervals(2)); %center
                
                fdomminb=fdommin(intervals(3):intervals(4)); %border
                fdomminm=fdommin(intervals(2):intervals(3)); %middle
                fdomminc=fdommin(intervals(1):intervals(2)); %center
               
                coveragedom=mean(f>=fdommin&f<=fdommax)*100;
                coveragedomc=mean(fc>=fdomminc&fc<fdommaxc)*100;
                coveragedomm=mean(fm>=fdomminm&fm<fdommaxm)*100;
                coveragedomb=mean(fb>=fdomminb&fb<fdommaxb)*100;
                
                %%%SUMMARY
                summary=[];
                summary.f=f;
                summary.fc=fc;
                summary.fm=fm;
                summary.fb=fb;
                summary.ffull=ffull;
                summary.ffullc=ffullc;
                summary.ffullm=ffullm;
                summary.ffullb=ffullb;
                
                %%%SUMMARY QUT
                summary.fmedian=fmedian;
                summary.fmax=fmax;
                summary.fmin=fmin;
                summary.fhat=fhat;
                
                summary.fmedianc=fmedianc;
                summary.fmaxc=fmaxc;
                summary.fminc=fminc;
                summary.fhatc=fhatc;
                
                summary.fmedianm=fmedianm;
                summary.fmaxm=fmaxm;
                summary.fminm=fminm;
                summary.fhatm=fhatm;
                
                summary.fmedianb=fmedianb;
                summary.fmaxb=fmaxb;
                summary.fminb=fminb;
                summary.fhatb=fhatb;
                
                summary.coverage=coverage;
                summary.coveragec=coveragec;
                summary.coveragem=coveragem;
                summary.coverageb=coverageb;
                
                summary.fmedianfull=fmedianfull;
                summary.fmaxfull=fmaxfull;
                summary.fminfull=fminfull;
                summary.fhatfull=fhatfull;
                
                summary.fmedianfullc=fmedianfullc;
                summary.fmaxfullc=fmaxfullc;
                summary.fminfullc=fminfullc;
                summary.fhatfullc=fhatfullc;
                
                summary.fmedianfullm=fmedianfullm;
                summary.fmaxfullm=fmaxfullm;
                summary.fminfullm=fminfullm;
                summary.fhatfullm=fhatfullm;
                
                summary.fmedianfullb=fmedianfullb;
                summary.fmaxfullb=fmaxfullb;
                summary.fminfullb=fminfullb;
                summary.fhatfullb=fhatfullb;
                
                summary.coveragefull=coveragefull;
                summary.coveragefullc=coveragefullc;
                summary.coveragefullm=coveragefullm;
                summary.coveragefullb=coveragefullb;
                
                
                %%%SUMMARY SA
                summary.fdommedian=fdommedian;
                summary.fdommax=fdommax;
                summary.fdommin=fdommin;
                summary.fhat=fhat;
                
                summary.fdommedianc=fdommedianc;
                summary.fdommaxc=fdommaxc;
                summary.fdomminc=fdomminc;
                summary.fhatc=fhatc;
                
                summary.fdommedianm=fdommedianm;
                summary.fdommaxm=fdommaxm;
                summary.fdomminm=fdomminm;
                summary.fhatm=fhatm;
                
                summary.fdommedianb=fdommedianb;
                summary.fdommaxb=fdommaxb;
                summary.fdomminb=fdomminb;
                summary.fhatb=fhatb;
                
                summary.coveragedom=coveragedom;
                summary.coveragedomc=coveragedomc;
                summary.coveragedomm=coveragedomm;
                summary.coveragedomb=coveragedomb;
                
                summary.fdommedianfull=fdommedianfull;
                summary.fdommaxfull=fdommaxfull;
                summary.fdomminfull=fdomminfull;
                summary.fhatfull=fhatfull;
                
                summary.fdommedianfullc=fdommedianfullc;
                summary.fdommaxfullc=fdommaxfullc;
                summary.fdomminfullc=fdomminfullc;
                summary.fhatfullc=fhatfullc;
                
                summary.fdommedianfullm=fdommedianfullm;
                summary.fdommaxfullm=fdommaxfullm;
                summary.fdomminfullm=fdomminfullm;
                summary.fhatfullm=fhatfullm;
                
                summary.fdommedianfullb=fdommedianfullb;
                summary.fdommaxfullb=fdommaxfullb;
                summary.fdomminfullb=fdomminfullb;
                summary.fhatfullb=fhatfullb;
                
                summary.coveragedomfull=coveragedomfull;
                summary.coveragedomfullc=coveragedomfullc;
                summary.coveragedomfullm=coveragedomfullm;
                summary.coveragedomfullb=coveragedomfullb;
                
                %%%%%
                
                %summary.msebootmean=sum((fdommean-f).^2/(M/2));
                %summary.msebootmedian=sum((fmedian-f).^2/(M/2));
                %summary.mse=sum((curresult.fhat'-f).^2/(M/2));
                %summary.logmsebootmean=sum((log(fmean)-log(f)).^2/(M/2));
                %summary.logmsebootmedian=sum((log(fmedian)-log(f)).^2/(M/2));
                %summary.logmse=sum((log(curresult.fhat)'-log(f)).^2/(M/2));
                save(filename,'curresult','sol','summary');
            end
        end
	end
end

