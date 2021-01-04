
load('data_and_results/realdata.mat')
M=512;
xqutfull=linspace(-M/2*sqrt(2),M/2*sqrt(2),M);
xqut=linspace(1,M/2*sqrt(2),M/2);
clear ftel

%Extract results from all realdata simulations
    for i=1:5
        i
        %Join all wavelets
        result=struct('Dhatfull',[],'Fhatfull',[],'Fhat',[],'Fhati',[],'Fhatd',[]);
        filelist=dir(strcat('data_and_results/*ind',num2str(i),'*.mat'));
        Ereal=realdata(i).E;
        Dhatfull_wavelet=zeros(length(filelist),M^2);
        
        for j=1:length(filelist)
            filename=filelist(j).name;
            load(filename)
            
            Fhatfull=sol.Fhatpos;
            Dhat=sol.Dhat;
            
            Dhatfull=zeros(size(Fhatfull,1),M^2);
            Fhat=zeros(size(Fhatfull,1),M/2);
            Fhati=zeros(size(Fhatfull,1),M/2);
            Fhatd=zeros(size(Fhatfull,1),M/2);
            
            if ~isempty(strfind(filename,'p'))
                Fhatfull=fliplr(Fhatfull);
            end

            %process results
            for k = 1:size(Fhatfull,1)
                fhat=Fhatfull(k,:);
                dhat=Dhat(k,:);
                dhat=reshape(dhat,528,528);
                dhat=dhat(1:512,1:512);
                
                if ~isempty(strfind(filename,'p'))
                    dhat=fliplr(dhat);
                end
                dhat=dhat.*Ereal(1:512,1:512);
                dhat=dhat(:)';

                fhati=fliplr(fhat(1:(M/2)));
				fhatd=fhat((M/2+1):end);
                fhatfull=fhat; 
                fhat=(fhati+fhatd)/2;
                
                fhati=interp1(xqut,fhati,1:(M/2),'linear');
                fhatd=interp1(xqut,fhatd,1:(M/2),'linear');
                fhat=interp1(xqut,fhat,1:(M/2),'linear');
                fhatfull=interp1(xqutfull,fhatfull,linspace(-M/2,M/2,M),'linear');
                
                Dhatfull(k,:)=dhat;
                Fhatfull(k,:)=fhatfull;
                Fhat(k,:)=fhat;
                Fhati(k,:)=fhati;
                Fhatd(k,:)=fhatd;
                
                if k==1 
                    Dhatfull_wavelet(j,:)=dhat;
                end
                    
            end
            
            result.Dhatfull_wavelet=Dhatfull_wavelet;
            result.Dhatfull=[result.Dhatfull;Dhatfull];
            result.Fhatfull=[result.Fhatfull;Fhatfull];
            result.Fhat=[result.Fhat;Fhat];
            result.Fhati=[result.Fhati;Fhati];
            result.Fhatd=[result.Fhatd;Fhatd];
            
        end
        
        result.fhat=median(result.Fhat);
        result.fhatd=median(result.Fhatd);
        result.fhati=median(result.Fhati);
        result.fhatfull=median(result.Fhatfull);
        
        Freal=realdata(i).F;
        Ereal=realdata(i).E;
        Oreal=realdata(i).O;
  
        %bootstrap for state of the art
        nboot=50;
        [Fdom,Fdomfull]=image_bootstrapSA(Freal,Ereal,Oreal,M,0.03,4,nboot);
        Fdomi=fliplr(Fdomfull(:,1:(M/2)));
        Fdomd=Fdomfull(:,(M/2+1):M);
        
        %Get estimation with state of the art
        solfull=astroStateOfTheArt(Freal,Ereal,Oreal, [0 360]);
        sold=astroStateOfTheArt(Freal,Ereal,Oreal, [0 180]);
        soli=astroStateOfTheArt(Freal,Ereal,Oreal, [180 360]);
        
        fdom=interp1(solfull.xbin,solfull.fbin,1:(M/2),'linear','extrap');
        fdom(fdom<=0)=min(solfull.fbin(solfull.fbin>0));
        
        fdomd=interp1(sold.xbin,sold.fbin,1:(M/2),'linear','extrap');
        fdomd(fdomd<=0)=min(sold.fbin(sold.fbin>0));
        
        fdomi=interp1(soli.xbin,soli.fbin,1:(M/2),'linear','extrap');
        fdomi(fdomi<=0)=min(soli.fbin(soli.fbin>0));
        
        fdomfull=[fliplr(fdomi) fdomd];
        
        result.fdom=fdom;
        result.fdomi=fdomi;
        result.fdomd=fdomd;
        result.fdomfull=fdomfull;
        result.Fdom=Fdom;
        result.Fdomfull=Fdomfull;
        result.Fdomi=Fdomi;
        result.Fdomd=Fdomd;
        result.fdomhat=median(Fdom);
        result.fdomhati=median(Fdomi);
        result.fdomhatd=median(Fdomd);
        result.fdomhatfull=median(Fdomfull);
        
        result.telescope=realdata(i).telescope;
        result.name=realdata(i).name; 
        
        ftel(i)=result;
    end
    
 
save('data_and_results/resultsRealData.mat','ftel');
    
 
    