function astroplot(sol,F,f,d)
    
    R=size(F,1)/2;
    %plot results
    ps=sol.options.ps;
    smoothf=sol.options.smoothf;
    dhat=sol.dhatpos;
    M=length(sol.fhat);
    n=size(F,1);
    
    
    %subplot(2^ps,2,1)
    
    imagesc(log(F'));
    
    title('Image')

    x=linspace(-M/2,M/2,M)*R/(M/2);
    x=x((M/2+1):M);
    
    logf=sol.logf;
    fhat=sol.fhat;
    
    logfi=logf((M/2):-1:1);
    logfd=logf((M/2+1):M);
    fhati=fhat((M/2):-1:1);
    fhatd=fhat((M/2+1):M);
    
    %subplot(2^ps,2,2:2:2^(ps+1))
    subplot(1,2,2)
    sol.cycle=0;

    %plot single
    curfhati=fhati;
    fhati(fhati<=0)=NaN;
    loglog(x,fhati,'color','blue','LineWidth',2)
    fhati=curfhati;
    hold on
    fhati(fhati>0)=NaN;
    %loglog(x,fhati,'*','color','blue')

    curfhatd=fhatd;
    fhatd(fhatd<=0)=NaN;
    loglog(x,fhatd,'color','red','LineWidth',2)
    
    %if(exist('sol.dominique')) 
%        loglog(sol.dominique.xbin,sol.dominique.fbin,'color','green','LineWidth',2); 
    %end
    
    fhatd=curfhatd;
    fhatd(fhatd>0)=NaN;
    %loglog(x,fhatd,'*','color','red')
    
    if sol.nboot~=0
    %plot cycles
        Fhat=sol.Fhat;
        
        Fhati=Fhat(:,(M/2):-1:1);
        Fhatd=Fhat(:,(M/2+1):M);
        
        fmean=mean([curfhatd;curfhati]);
        Fmax=fmean-quantile([Fhati;Fhatd],0.95);
        Fmin=fmean-quantile([Fhati;Fhatd],0.05);
        
        errorbar(x,mean([curfhatd;curfhati]),Fmin,Fmax,' ','color',[0.5 0.5 0.5]);
        loglog(x,curfhati,'color','blue','LineWidth',2)
        loglog(x,curfhatd,'color','red','LineWidth',2)
        
%         curFhati=Fhati;
%         Fhati(Fhati<=0)=NaN;
%         loglog(x,Fhati',':','color', 0.7*[1 1 1])
%         
%         
%         curFhatd=Fhatd;
%         Fhatd(Fhatd<=0)=NaN;
%         loglog(x,Fhatd',':','color', 0.7*[1 1 1])
        
    end
%    loglog(sol.dominique.xbin,sol.dominique.fbin,'color','green','LineWidth',2); 
    %plot true profile if known
    if nargin>2 
       
        xf=linspace(-n/2,n/2,n)*R/(n/2);
        xf=xf((n/2+1):n);
        fi=f((length(f)/2):-1:(length(f)/2-n/2+1));
        fd=f((length(f)/2+1):(length(f)/2+n/2));
        loglog(xf,fi,'color',[0 0 0]);
        if sum(fi~=fd)~=0
            loglog(xf,fd,'color',[0.5 0.5 0.5]);
        end
    end
    hold off
    title('Profile function')

    %plot point sources
    if ps
        %subplot(2,2,3)
        subplot(1,2,1)
        imagesc(log(F'));
        title('Detected Point Sources')
        hold on

        if nargin>3
            %plot true point sources if known
            indd=find(d~=0);
            xd=mod(indd-1,n)+1;
            yd=ceil(indd/n);
            plot(xd,yd,'ro','MarkerSize',6)
        end

        %if sol.nboot==0
            indd=find(dhat~=0);
            xd=mod(indd-1,n)+1;
            yd=ceil(indd/n);
            plot(xd,yd,'wo','MarkerSize',6)

             if nargin>3
                %plot true point sources if known
                indd=find(d'==0 & dhat~=0);
                xd=mod(indd-1,n)+1;
                yd=ceil(indd/n);
                plot(xd,yd,'w*','MarkerSize',4)
             end
        %end
        
        
        hold off    
    end
    %qmf=MakeONFilter('Daubechies',8);
    %fw=FWT_PO(f',0,qmf);
    %figure
    %plot(fw,'.k' )
    %hold on;
    %plot(sol.fwlasso,'.b');
    %plot(sol.fwhat, '.r');
    %hold off
end