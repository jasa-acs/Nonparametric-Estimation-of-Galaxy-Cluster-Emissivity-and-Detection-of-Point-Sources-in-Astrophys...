%function [sol, F, f, d] = astrosimPS(k,n,R,sf,type,blurtype,wavelets,nitr, O,E,ps,psize,pnum,rotate,smoothf,showplot,parallel)
function [sol, F, f, d] = astrosimPS(k,n,R,E,O,sf,type,psize,nboot,parallel,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulate image for a given M and find profile function and point sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%k       -> Seed for random generation. Set 0 for completely random (no seed)  
%M       ->	Generate a MxM image
%R       ->	Radius of the cluster (length between the center to the border of
%           the image
%O       ->	Background. It must be in vector form of size of a MxM. 
%E       ->	Sensitivity. 
%sf      -> Scale factor for the profile
%type    -> profile shape. Options are 'rho', 'ex', 'null'.
%psize  -> Size of the simulated point sources
%nboot   -> Number of bootstraps (set 0 for no bootstrapping);
%parallel-> 1 for doing in parallel

%options -> SEE  help astrosolvePS"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
%solution -> structure with the following elelments:
%
%logf     -> Estimated log profile
%fhat     -> Estimated profile
%fwhat    -> Estimated wavelet coefficients
%dhat     -> Estimated point sources
%F        -> Simulated Image
%f        -> Simulated profile
%d        -> Simulated point sources

%If rotate=1 then
%Fhat     -> Estimated profile of all wavelets and rotations
%Fwhat    -> Estimated wavelet coefficients using GLM of all wavelets and rotations
%Dhat     -> Estimated point sources for all wavelets and rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pnum=100;
    %Simulate profile
    if k==0 
        k=floor(rand(1)*1000)+1; 
    end
    rng(k); %seed
    
    f=generatef(n,sf,type);
    N=n^2;
    %Create Abel Matrix
    [K,B,qmf,S]=setUpOperatorsWS(n,N,R,options.bluralpha,options.blurR0,options.blurthresh,[],[],1,1);

    %simulate point sources
    d=zeros(1,N);
    
    %simulate image
    mu0=K*f';
    if options.ps
        %Fix point sources
        indd=floor([n/2+n/4+1 n/2; n/2-n/4 n/2+n/4+1;   n/2+3*n/8+1 n/2-n/4;         n/2-3*n/8+1 n/2-n/8;         n/2-n/8+1 n/2+n/8+1]);
        d=zeros(n,n);
        for i=1:5
            d(indd(i,1),indd(i,2))=psize;
        end
        d=d(:)';
        
        %random point sources
        imind=reshape(1:N,sqrt(N),sqrt(N));
        nops=imind((floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)),(floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)));
        clear imind
        nops=nops(:);
        yesps=setdiff(1:N,nops);
        d=d*0;
        numofps=sqrt(N)/4;
       
       

        %POWER LAW
        randalpha=1.82;
        %psize=max(mu0);
        Smax=mean(mu0)*500;
        %Smin=mean(mu0)*50;
         
        Cpowerlaw=1/(Smax^-randalpha);
        dminpowerlaw=(N/Cpowerlaw)^(-1/randalpha);
        Smin=2*dminpowerlaw;
        
        dsizes=[0 linspace(Smin,Smax,N)];
        %dsizes=[0 Smin:(2*dminpowerlaw):Smax];
        Fpowerlaw=Cpowerlaw*dsizes.^-randalpha;
        Fpowerlaw(1)=N;
        Fpowerlaw=[Fpowerlaw,0];
        fpowerlaw=Fpowerlaw(1:(length(Fpowerlaw)-1))-Fpowerlaw(2:length(Fpowerlaw));
        d = dsizes(randsample(length(dsizes),N,true,fpowerlaw));
        [dy,di]=sort(d,'descend');
        d=d*0;
        numofps=sqrt(N)/2;
        d(di(1:numofps))=dy(1:numofps);
        
%         curDnum=0;
%         numofps=N
%         for ii=length(dsizes):-1:1
%             dnum=floor(Cpowerlaw*dsizes(ii)^-randalpha-curDnum);
%             if (dnum>0 & curDnum<numofps)
%                 ii
%                 floor(Cpowerlaw*dsizes(ii)^-randalpha)
%                 dnum=min(dnum,numofps-curDnum);
%                 whichd=randsample(find(d==0),dnum);
%                 d(whichd)=dsizes(ii);
%                 curDnum=curDnum+dnum;
%             end
%         end
        
        %randk=rand(1,length(d));
        %power_law=(1-randk).^(-1/randalpha)-1;
        %d=power_law*5e-7;
        
        %d=d*0;d(randsample(yesps,floor(numofps)))=rand(1,numofps)*psize;
    end
    nonoise=0;
    
    d(mu0==0)=0;
    mu0=mu0+d';
    if length(E)==1
        if(E==0)
            nonoise=1;
            E=1;
        end
        E=E*ones(n,n);
    else
        E=imresize(E,[n n],'method','nearest');
    end
    E=E(:);
    E(mu0==0)=0;
    
    %add blurring
    mu0=reshape(mu0',n,n);
    mu0=filter2(B,mu0);
    mu0=mu0(:);
    
    %add sensitivity and offset
    if length(O)==1
        O=O*ones(n,n);
    end
    O=O(:);
    O=O-min(min(O))+1e-4;
    O(mu0==0)=1e-4;
    
    mu0=E.*mu0+O;
    
    %Add poisson noise
    F=0;
    safecount=0;
    while (sum(F)==0&safecount<100)
        safecount=safecount+1;
        'Adding Poisson noise'
        if nonoise==0
            F = poissrnd(mu0',1,N);
        else
            F=mu0';
        end
    end
    
    F=reshape(F,n,n);
    E=reshape(E,n,n);
    O=reshape(O,n,n);

%    save(strcat('EOsim',num2str(n),'.mat'),'E','O')
%imagesc(log(F));
%    asdasdasd
    %get dominique's current results

    %dominique=astrodominique(F,E,O);
    dominique=0;
    
    %sol=F;f=f;d=0;return
    showplot=options.showplot;
    options.showplot=0; %Don't plot single process but the full simulation
    if nboot==0
        %[sol,F]= astrosolvePS(F,R,E,O,[],options);
        %[sol,F]= astrosolvespline(F,R,E,O,[],options);
        [sol,F]= astrosolveWS(F,R,E,O,[],options);
    else
        %[sol,F]= astrobootPS(F,R,E,O,nboot,options);
        [sol,F]= astrobootWS(F,R,E,O,nboot,options);
    end
    d=reshape(d,n,n);
    newn=sqrt(length(sol.dhat));
    d=d((n/2-newn/2+1):(n/2+newn/2),:);
    d=d(:,(n/2-newn/2+1):(n/2+newn/2));
    d=d(:)';

    sol.nboot=nboot;
    sol.dominique=dominique;
    %plot
    if (~parallel & showplot)
        if(strcmp(type,'blocks'))
            plotblock
        else
            astroplot(sol,F,f,d);
        end
    end

    sol.R=R;
end
