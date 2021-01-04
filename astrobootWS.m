function [sol,F]= astrobootWS(F,R,E,O,nboot,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%F       ->	Image of size nxn to estimate profile and point sources.
%R       ->	Radius of the cluster (length between the center to the border of
%           the image
%nboot   -> Number of bootstraps samples
%O       ->	Background. It must be in image form of size of nxn. If O=0 then 
%           there is no offset. If O is a constant, then is constant offset.
%E       ->	Sensitivity. It must be in image form of size of a nxn. If E=1 then 
%           there is no sensitivity effect. If E is a constant, then is constant sensitivity.
%Options ->


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
%solution -> structure with the following elelments:
%
%logf     -> Estimated log profile
%fhat     -> Estimated profile
%fhatpos  -> Estimated profile with positive values
%fwhat    -> Estimated wavelet coefficients
%dhat     -> Estimated point sources
%Fhat     -> Estimated profile of bootstraps
%Fhatpos  -> Estimated profile of bootstraps with just positive values
%Fwhat    -> Estimated wavelet coefficients for all bootstraps
%Dhat     -> Estimated point sources for all bootstraps
%Flasso   -> Estimated profile of bootstraps with lasso coefficients
%Fwlasso  -> Estimated wavelet lasso coefficients for all bootstraps
%Dlasso   -> Estimated lasso point sources for all bootstraps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parallel=options.parallel;
    %Check sizes, make square, etc.
    [F,E,O]=astro_setFEO(F,E,O,0,options.remblack);
    oldn=ceil(size(F,1));
    n=oldn-ceil(options.center*oldn);
    M=2^floor((log(n)/log(2)));
    
    waveletoptions=[2 4; 2 6; 2 8; 1 4; 1 6; 1 8; 3 2; 3 3; 3 4];
    
    %masksize=ceil(sqrt((options.blurthresh^(-1/options.bluralpha)-1)*options.blurR0^2));
    %masksize=2*masksize+1;
    %bootwindow=(2^floor(log(masksize)/log(2)))
    
    masksize=ceil(sqrt((options.blurthresh^(-1/options.bluralpha)-1)*options.blurR0^2));
    ii=[masksize:-1:0 1:masksize];
    xx=repmat(ii,length(ii),1);
    yy=repmat(ii',1,length(ii));
    d=sqrt(xx.^2+yy.^2);
    %bootwindow=(1+(d/options.blurR0).^2).^(-options.bluralpha);
    %bootwindow=bootwindow/sum(sum(bootwindow));
    
    MM=0;
    if(options.sshrink>0) MM=MM+M; end
    if(options.wavelet(1)>=0) MM=MM+M; end

    %Create return vectors
    Logf=zeros(nboot,M);
    Fhat=zeros(nboot,M);
    Fhatpos=zeros(nboot,M);
    FWhat=zeros(nboot,MM);
    Dhat=zeros(nboot,oldn^2);
    Flasso=zeros(nboot,M);
    FWlasso=zeros(nboot,MM);
    Dlasso=zeros(nboot,oldn^2);

    K=abelmatrix2Half2D(M,n^2,R,1,1,1);   %revisar
    options.showplot=0;
    options.remblack=0;
    
    
    if(~parallel)
        for i=1:nboot
            strcat('nboot',num2str(i))
            options1=options;
            waveletoptions1=waveletoptions;
            
            %bootstrap image in 2x2 squares
            if i==1 
                %first do pure image (no bootstrap)
                F1=F(((oldn-M)/2+1):((oldn-M)/2+M),((oldn-M)/2+1):((oldn-M)/2+M));
                E1=E(((oldn-M)/2+1):((oldn-M)/2+M),((oldn-M)/2+1):((oldn-M)/2+M));
                O1=O(((oldn-M)/2+1):((oldn-M)/2+M),((oldn-M)/2+1):((oldn-M)/2+M));
                up_move=0; right_move=0;
            else
                %boottrap image
                [F1,E1,O1, up_move,right_move]=image_bootstrap(F,E,O,options.center,options1.bootwindow);
            end
            
            %choose wavelet
            if(options1.wavelet(1)~=-1)
                if length(options1.wavelet)==1
                    options1.wavelet=waveletoptions1(options1.wavelet,:);
                else
                    options1.wavelet=waveletoptions1(mod(i,9)+1,:);
                end
            end
            
            %solve single bootstrap
            sol1= astrosolveWS(F1,R,E1,O1,[],options1);

            Logf(i,:)=sol1.logf;
            Fhat(i,:)=sol1.fhat;
            Fhatpos(i,:)=sol1.fhatpos;
            FWhat(i,:)=sol1.fwhat;
            Flasso(i,:)=sol1.flasso;
            FWlasso(i,:)=sol1.fwlasso;
            
            if(options.center>0)
                %restore position of point sources
                Dlasso(i,:)=rewindMoveCenter(sol1.dlasso,n,size(F,1),up_move,right_move);
                Dhat(i,:)=rewindMoveCenter(sol1.dhat,n,size(F,1),up_move,right_move); 
            else
                Dhat(i,:)=sol1.dhat;
                Dlasso(i,:)=sol1.dlasso;
            end
        end
    elseif(parallel)
         parfor i=1:nboot
            
            strcat('nboot',num2str(i))
            options1=options;
            waveletoptions1=waveletoptions;
            
            %bootstrap image in (bootwindow x bootwindow) squares and changing center (if options.center>0)
            if i==1 
                %first do pure image (no bootstrap)
                F1=F;E1=E;O1=O;
                F1=F1(((oldn-M)/2+1):((oldn-M)/2+M),((oldn-M)/2+1):((oldn-M)/2+M));
                E1=E1(((oldn-M)/2+1):((oldn-M)/2+M),((oldn-M)/2+1):((oldn-M)/2+M));
                O1=O1(((oldn-M)/2+1):((oldn-M)/2+M),((oldn-M)/2+1):((oldn-M)/2+M));
                up_move=0; right_move=0;
            else
                %boottrap image
                [F1,E1,O1, up_move,right_move]=image_bootstrap(F,E,O,options.center,options1.bootwindow);
            end
            
            %choose wavelet
            if length(options1.wavelet)==1
                options1.wavelet=waveletoptions1(options1.wavelet,:);
            else
                options1.wavelet=waveletoptions1(mod(i,9)+1,:);
            end
            options1
            
            %solve single bootstrap
            sol1= astrosolveWS(F1,R,E1,O1,[],options1);

            Logf(i,:)=sol1.logf;
            Fhat(i,:)=sol1.fhat;
            Fhatpos(i,:)=sol1.fhatpos;
            FWhat(i,:)=sol1.fwhat;
            Flasso(i,:)=sol1.flasso;
            FWlasso(i,:)=sol1.fwlasso;
           
            if(options.center>0)
                %restore position of point sources
                Dlasso(i,:)=rewindMoveCenter(sol1.dlasso,n,size(F,1),up_move,right_move);
                Dhat(i,:)=rewindMoveCenter(sol1.dhat,n,size(F,1),up_move,right_move); 
            else
                Dhat(i,:)=sol1.dhat;
                Dlasso(i,:)=sol1.dlasso;
            end

        end    
    end
    
    clear E
    clear O
    
    %Clean connected point sources
    if nboot~=1
        dpos=max(Dhat); 
    else
        dpos=Dhat;
    end
    
    if  options.ps
        %choose center of connected point sources
        dpos=reshape(dpos,oldn,oldn);
        CC = bwconncomp(dpos>0);
        for i=1:CC.NumObjects
            dpos(CC.PixelIdxList{i})=dpos(CC.PixelIdxList{i}).*(dpos(CC.PixelIdxList{i})==max(dpos(CC.PixelIdxList{i})));
        end
        dpos=dpos(:)';
     
        %filter small values      
        dpostemp=dpos/max(dpos);
        dpos=dpos.*(dpostemp>min(graythresh(dpostemp),max(dpostemp)/100));
        
        %extrapolate
        dpos=reshape(dpos,oldn,oldn);
        dpos=imresize(dpos,size(F),'nearest');
        CC = bwconncomp(dpos>0);
        for i=1:CC.NumObjects
            dpos(CC.PixelIdxList{i})=dpos(CC.PixelIdxList{i}).*(F(CC.PixelIdxList{i})==max(F(CC.PixelIdxList{i})));
        end
        dpos=dpos(:)';
        
    end
    
    if nboot~=1
        fhat=median(Fhat);
        fhatpos=median(Fhatpos);
    else
        fhat=Fhat;
        fhatpos=Fhatpos;
    end
    if nboot~=1
        flasso=median(Flasso);
    else
        flasso=Flasso;
    end
    sol.flasso=flasso;
    sol.Flasso=Flasso;
    sol.fhat=fhat;
    sol.fhatpos=fhatpos;
    sol.fwhat=median(FWhat);
    sol.dhat=dpos';
    sol.Logf=Logf;
    sol.Fhat=Fhat;
    sol.Fhatpos=Fhatpos;
    sol.FWhat=FWhat;
    sol.Dhat=Dhat;
    sol.nboot=nboot;
    sol.dhatpos=sol.dhat;
    
    %get log(f) and interpolate negative values of f
    y=log(fhat);
    y(fhat<0)=NaN;
    x=linspace(-M/2,M/2,M);
    xi=x(find(~isnan(y)));yi=y(find(~isnan(y)));
    logf=interp1(xi,yi,x,'linear','extrap');

    %smooth result
    if options.smoothf  logf=smooth(logf); end;
    sol.logf=logf;
    
    sol.options=options;
    if options.showplot astroplot(sol,F); end;
    sol.R=R;
end 

function [F,E,O,up_move,right_move]=image_bootstrap(F,E,O,center,bootwindow)
    n=size(F,1);
    if center>0
        newn=n-ceil(center*n);
        up_move=floor(rand(1)*ceil(center*n));
        right_move=floor(rand(1)*ceil(center*n));
        F=F((1+up_move):(up_move+newn),(1+right_move):(right_move+newn));
        E=E((1+up_move):(up_move+newn),(1+right_move):(right_move+newn));
        O=O((1+up_move):(up_move+newn),(1+right_move):(right_move+newn));
    else
        up_move=0;right_move=0;
    end
    n=size(F,1);
 
%% 2 points symmetric bootstrap
    for i=1:(n/2)
        for j=1:n
            symInd=[i,(n-i+1)];
            bootInd=randsample(symInd,2,true);
            F(symInd,j)=F(bootInd,j);
            E(symInd,j)=E(bootInd,j);
            O(symInd,j)=O(bootInd,j);
        end
    end
    

    
%% blurring mask window

%     F1=F;E1=E;O1=O; 
%     masksize=(size(bootwindow,1)-1)/2;
% 
%     for i=(masksize+1):(n-masksize)
%         for j=(masksize+1):(n-masksize)
%             i0=i-masksize;
%             i1=i+masksize;
%             j0=j-masksize;
%             j1=j+masksize;
%             
%             Fwindow=F1(i0:i1,j0:j1);
%             Ewindow=E1(i0:i1,j0:j1);
%             Owindow=O1(i0:i1,j0:j1);
% 
%             ijwindow=randsample(1:size(bootwindow,1)^2,1,true,bootwindow(:));
%             
%             F(i,j)=Fwindow(ijwindow);
%             E(i,j)=Ewindow(ijwindow);
%             O(i,j)=Owindow(ijwindow);
% 
%         end
%     end

%% nxn window bootstrap

    if(bootwindow==0)
        bootwindow=floor(n/64);
        strcat('bootwindow',num2str(bootwindow))
    end
    
    for i=1:(n/bootwindow)
        for j=1:(n/bootwindow)
            i0=bootwindow*(i-1)+1;
            i1=bootwindow*i;
            j0=bootwindow*(j-1)+1;
            j1=bootwindow*j;
            
            Fwindow=F(i0:i1,j0:j1);
            Ewindow=E(i0:i1,j0:j1);
            Owindow=O(i0:i1,j0:j1);
            ijwindow=randsample(1:bootwindow^2,bootwindow^2,true);
            
            Fwindow(:)=Fwindow(ijwindow);
            Ewindow(:)=Ewindow(ijwindow);
            Owindow(:)=Owindow(ijwindow);
            
            F(i0:i1,j0:j1)=Fwindow;
            E(i0:i1,j0:j1)=Ewindow;
            O(i0:i1,j0:j1)=Owindow;
        end
    end
    
end

function D=rewindMoveCenter(D,newn,oldn,up_move,right_move)
    D=reshape(D,newn,newn);
    Dtemp=zeros(oldn,oldn);
    Dtemp((1+up_move):(up_move+newn),(1+right_move):(right_move+newn))=D;
    D=Dtemp(:)';
end



    
    
    
    
        
    