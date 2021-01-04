function [sol,F]= astrosolveWS(F,R,E,O,K,options)
%
%Finds the profile function and point sources in an image of size n x n. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL
%
%                 F~Pois(E.B*[KW I][beta d]+O)
%
%where K is the Abel Transform, B is the blurring, E is the sensitivity, O is the offset 
%f=W beta is the profile and d is the point source matrix. W is wavelet
%matrix, beta is wavelet coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%F       ->	Image to estimate profile and point sources, must be square image!
%R       ->	Radius of the cluster (length between the center to the border of the image
%E       ->	Sensitivity. It must be in image form of size of a nxn. If E=1
%then  
%           there is no sensitivity effect. If E is a constant, then is constant sensitivity.
%O       ->	Background. It must be in image form of size of nxn. If O=0 then 
%           there is no offset. If O is a constant, then is constant offset.
%K       -> K matrix (computed if length(k)==1)

%options ->Structure with the following options. Check defaults with astro_setparams()

    %bluralpha->  parameter for psf (Set to 0 for no blurring)  psf=point spread function
    %blurR0-> parameter for psf
    %blurthresh->  parameter for psf
    
    %wavelet -> 2 entries vector. First entry is the wavelet type and second is
    %           wavelet parameter.  Wavelet type is 1 for 'Daubechies', 
    %           2 for 'Symmlet', 3 for 'Coiflet'. If length(wavelet)=1 then
    %           wavelet is set to [1 8].
    %nitr    -> Number of iterations for convergence, nitr=200 is enough but
    %           nitr=1000 is better to be sure
    %ps      -> If there are point sources in image, set to 1.  Set to 0 if not
    %psclean -> Set to 1 if point sources cleaning is desired after FISTA and GLM
    %smoothf -> If a smoothing filter should be applied after profile estimation
    %showplot-> Show plot
    %prevfw  -> value of not penalized wavelet coefficients (taken into account
    %           if length(prevfw)~=0))
	%MORE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
%Structure with the following elemets:
%
%logf     -> Estimated log profile
%fhat     -> Estimated profile
%fhatpos  -> Estimated profile with forced positive values
%fwhat    -> Estimated wavelet coefficients using GLM
%dhat     -> point sources using GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options0=options;
    % Check number of arguments
    if(nargin < 5 || nargin> 6 )
    error('incorrect number of arguments')
    end
    
    % Set default parameters or grab user-defined params
    if (nargin < 6 || isempty(options) )
       opts = astro_setparams();
    else
       opts = astro_setparams(options);
    end
    
    blurthresh =opts.blurthresh;
    bluralpha =opts.bluralpha;
    blurR0    =opts.blurR0;
    wavelet   =opts.wavelet;
    nitr      =opts.nitr;
    ps        =opts.ps;
    psclean   =opts.psclean;
    smoothf   =opts.smoothf;
    showplot  =opts.showplot;
    prevfw    =opts.prevfw;
    remblack  =opts.remblack;
    sshrink   =opts.sshrink;
    beta0     =opts.beta0;
    itersol   =opts.itersol;
    forcepos  =opts.forcepos;
    prevPS    =opts.prevPS;

    %Check Image size, make them squares of even size, convert to vectors
    [F,E,O]=astro_setFEO(F,E,O,1,opts.remblack);
    R=size(F,1)/2;
    
    %dominique=astrodominique(F,E,O);
    %Convert images in vectors
    %F(F==0)=1e-10;
%     subplot(1,2,1) 
%     imagesc(F)
%     subplot(1,2,2) 
%     imagesc(F>0)
%     sddsf
    F=F(:)';
    O=O(:)';
    E=E(:)';
    
    N=length(F);
    n=sqrt(N);
    M=2^floor(log(sqrt(N))/log(2));

    Lf=[];
    
    %weights, not used
    %w=ones(1,length(F))';
    
    %set wavelet type and parameter
    
    wavetype=wavelet(1);
    wavepar= wavelet(2);
    if(wavelet(1)<0)
        wavetype='null';
    elseif wavelet(1)==0
       wavetype='Haar'; 
    elseif wavelet(1)==1
        wavetype='Daubechies';
    elseif wavelet(1)==2
        wavetype='Symmlet';
    elseif wavetype==3
        wavetype='Coiflet';
    end
    %if(length(wavelet)<2) wavelet=[1 8]; end
    
    

    %Set up operators
    if isempty(K)
        [K,B,qmf,S]=setUpOperatorsWS(M,N,R,bluralpha,blurR0,blurthresh,wavetype,wavepar,sshrink,1);
    else
        [K,B,qmf,S]=setUpOperatorsWS(M,N,R,bluralpha,blurR0,blurthresh,wavetype,wavepar,sshrink,K);
    end
    
    MM=M*(~isempty(qmf))+M*(~isempty(S));
    
    % FISTA OPTIONS
    options.fistaitr  = nitr;
    options.fistatol  = 1e-6;
    options.fistalogperiod   = 1;
    
    if(forcepos==2) %Force all Kings
        options.indpos=(MM-M*(~isempty(S))+1):(MM+N*ps);
        if(~isempty(options.indpos))
            if(options.indpos(1)~=1) options.indpos=[1 options.indpos]; end;
        else
            options.indpos=1;
        end
    elseif(forcepos==1) %Force intercep
        options.indpos=1;
    else
        options.indpos=[]; %no force
    end

        
    options.ind1    = 1;
    %options.ind0    = find(opts.fistax0==0);
    
    %Estimate intercept
    I=zeros(MM,1);    I(1)=1;

    
    
    A=abeloperatorWS(I,1,K,B,qmf,S,ones(1,MM),E,0); 
    clear I
    if isempty(beta0)
        %beta0 is the estimate of the father wavelet coefficient under the null model
        if ps
            beta0=glmfit(A,F','poisson','constant','off','link','identity','offset',O) 
            beta0=sum(F~=0)*median(F(F~=0))/(sum(A))
        else
            beta0=glmfit(A,F','poisson','constant','off','link','identity','offset',O)
        end
        clear A
        itersolver=0;  
    else
        if (itersol~=0)
            itersolver=1;
        else
            itersolver=0;
        end
        
    end
    
    %beta0=1e-40
    %beta0=0
    %beta0=beta0/100
%         if(isempty(qmf))
%             beta0=0
%         end

    
    clear I
    clear A
    
    x0initial=prevfw;
    prevfw=0*ones(MM+N*ps,1);
    prevfw(1)=beta0;
    %Get Lambda and rescaling factors
    indx0initial=find(x0initial==0);
    if(beta0<0)
        save('badF.mat','F')
        beta0=0;
        prevfw(1)=beta0;
        x0initial=prevfw;
        warning('Negative intercept, setting to zero')
    end
    [lambda, divX,lambdaMax,divX1,prevPS]=lambdaQUTastroWS(F,K,prevfw,O,qmf,S,B,E,ps,prevPS);
    
    
    if (itersolver==0) 
        x0initial=prevfw;
        %warning('CUIDADO, ELIMINAR ESTO DESPUÉS')
    else
        divX(indx0initial)=Inf;
    end
    %SOLVE GLM L1 
    
    %Atemp=zeros(N,M);
    %for i = 1:M
    %    I=zeros(M,1);I(i)=1;
    %    Atemp(:,i)=abeloperatorPS(I,1,K,B,qmf,ones(1,M),E,0); 
    %end
    %betaini=glmfit(Atemp,F','normal','constant','off','link','identity','offset',O,'weights',w');

    %load beta0.mat
    %prevfw=beta0;
    %prevfw(1)=0

    options.fistax0=prevfw.*divX';
    options.fistax0(isnan(options.fistax0))=0;

    A  = operator(@(x,mode) abeloperatorWS(x,mode,K,B,qmf,S,divX,E,ps),N,MM+N*ps);
    A1=A;
    optionsF=options;
    sol.static  = FISTA_Poi(A, O', F', lambda, Lf , options,qmf,S,divX,lambdaMax);
    if(isempty(sol.static.x)) sol.static.x=options.fistax0; end
    sol.static.x=sol.static.x./divX';

    %setup estimates
    sol.glm=sol.static.x*0; 
    fwlasso=sol.static.x(1:MM);
    sol.static.x(1:MM)=[];
    
    dlassopure=sol.static.x;
    dlasso2=sol.static.x;
    dlasso3=dlasso2;

    %Detect point sources
    if ps
        dlasso2(dlasso2<0)=0;
        %dlasso3=reshape(dlasso2',n,n);
        %dlasso3=imdilate(dlasso2>0,ones(3,3));
        %dlasso3=dlasso3(:);
        if psclean
            dlasso2=reshape(dlasso2',n,n);
            CC = bwconncomp(dlasso2>0);
            for i=1:CC.NumObjects
                dlasso2(CC.PixelIdxList{i})=dlasso2(CC.PixelIdxList{i}).*(dlasso2(CC.PixelIdxList{i})==max(dlasso2(CC.PixelIdxList{i})));
            end
            dlasso2=dlasso2(:);
        end
    end
    
    sol.static.x= [fwlasso;dlasso2];
    
    %nnz=sum(sol.static.x~=0);
    %indnz=find(sol.static.x~=0);
    indnz=find([fwlasso; dlasso2]~=0);
    
    clear dlasso2
    nnz=length(indnz);
    
    %GET GLM FIT
        %build A matrix
        A=zeros(N,nnz);
        for i = 1:nnz
            I=zeros(MM+N*ps,1);
            I(indnz(i))=1;
            A(:,i)=abeloperatorWS(I,1,K,B,qmf,S,divX,E,ps);
        end
        
        %get beta from glm with MATLAB
        %betanew=glmfit(A,F','poisson','constant','off','link','identity','offset',O','weights',w);
        
        %if negative objective function, get with FISTA
        %if (imag(sum(A*betanew+O' - F' .* log(A*betanew+O')))~=0|1==1)
        	%options.fistaitr  = 1000;
            
            options.fistax0=sol.static.x(indnz).*divX(indnz)';
            indnotzero=find(sol.static.x~=0);
            %options.fistax0=sol.static.x(indnz)*0;
            %options.fistax0(1)=sol.static.x(1)*divX(1);
            options.indpos=[];
            if(forcepos==2)
                options.indpos=find(indnotzero>=(MM-M*(~isempty(S))+1));
                if(~isempty(options.indpos))
                    if(options.indpos(1)~=1) options.indpos=[1; options.indpos]; end;
                else
                    options.indpos=1;
                end
            elseif(forcepos==1)
                options.indpos=1;
            else
                options.indpos=[];
            end
            sol.glm0  = FISTA_Poi(A, O', F', 0, Lf , options,qmf,S,divX(indnz),0);
            if sol.glm0.conv==0
            	betanew=sol.static.x(indnz);
            else
                betanew=sol.glm0.x;   
            end
        %end    
        
    
    %get estimates
    
    sol.glm(indnz)=betanew./divX(indnz)';
    fwhat=sol.glm(1:MM);
    sol.glm(1:MM)=[];
    dhat=sol.glm;
    if ps==0 dhat=zeros(1,N); end
    
    fhats=0;flassos=0;fhatw=0;flassow=0;
    
    fwhattemp=fwhat;
    fwlassotemp=fwlasso;
    if(~isempty(qmf))
        fhatw=IWT_PO(fwhattemp(1:M),0,qmf);
        flassow=IWT_PO(fwlassotemp(1:M),0,qmf);
        fwhattemp(1:M)=[];
        fwlassotemp(1:M)=[];
    end
    if(~isempty(S))
        fhats=S*fwhattemp;
        flassos=S*fwlassotemp;
    end
    
    fhat=fhats+fhatw;
    flasso=flassos+flassow;

    sol.fhat=fhat;
    sol.fwhat=fwhat;

    sol.fwlasso=fwlasso;
    sol.flasso=flasso;
    if ps==0 dlassopure=dhat; end
    sol.dlasso=dlassopure;
    sol.dhat=dhat;
    sol.dhatpos=dhat*0;

    %interpolation-extrapolation
    
    logf=fhat*0+1;
    sol.logf=logf;
    fhatpos=fhat;
    if(~isempty(S))
        fhatpos(fhat<=0)=fhats(fhat<=0);
        fhatpos(fhatpos<=0)=min(fhatpos(fhatpos>0));
    else
        xx=linspace(-M/2,M/2,M);
        fhatpos(fhat<0)=NaN;
        logf=log(fhatpos);
        if(sum(fhat>0)<4)
            if(sum(fhat>0)<1) 
                fhatpos=log(fhat*0+1e-60);
            else
                logf(fhat<0)=min(fhat(fhat>0));
                fhatpos=exp(logf);
            end
        else
            fhatpos = exp(interp1(xx(fhat>0),logf(fhat>0),xx,'spline','extrap'))';
        end
        fhatpos(fhatpos<=0)=min(fhatpos(fhatpos>0));
    end
    
    sol.fhatpos=fhatpos;
    sol.options=options;

    if showplot astroplot(sol,reshape(F,sqrt(N),sqrt(N))); end;
    sol.R=R;
    sol.nboot=0;
    F=reshape(F,sqrt(length(F)),sqrt(length(F)));
    
    [beta0 sol.fwhat(1)]
    newbeta0=sol.fwhat(1);
    if(newbeta0<0) newbeta0=0; end;
    beta0error=abs(beta0-newbeta0)/beta0
    if(newbeta0==0) beta0error=0; end;
    
    if(beta0error>0.1&beta0~=newbeta0&((itersol<2&~ps)|(itersol<2&ps))&1==1)
        fhat0=sol.fhat;
        flasso0=sol.flasso;
        fwhat0=sol.fwhat;
        fwlasso0=sol.fwlasso;
        options=options0;
        options.beta0=sol.fwhat(1);
        options.prevfw=sol.fwlasso;
        options.prevfw(1:M)=options.prevfw(1:M)+0.00001
        if ps
            options.prevfw=[sol.fwlasso; sol.dlasso];
        end
        options.itersol=itersol+1;
        options.prevPS=prevPS;
        [sol,F]= astrosolveWS(reshape(F,n,n),R,reshape(E,n,n),reshape(O,n,n),K,options);
        sol.flasso0=flasso0;
        sol.fwhat0=fwhat0;
        sol.fwlasso0=fwlasso0;
    end
    
    %Check point sources
    if(ps==1&itersol==0)
    	O=abeloperatorWS(sol.fwhat,1,K,B,qmf,S,ones(1,MM),E,0)'+O; 
        O(O<=0)=0.0001;
        dhattemp=reshape(sol.dhat',n,n);
        CC = bwconncomp(dhattemp>0);
        for i=1:CC.NumObjects
            dhattemp(CC.PixelIdxList{i})=dhattemp(CC.PixelIdxList{i}).*(dhattemp(CC.PixelIdxList{i})==max(dhattemp(CC.PixelIdxList{i})));
        end
        dhattemp=dhattemp(:);
        
%         indnz=MM+find(dhattemp~=0);
%         nnz=length(indnz);
% 
%         %GET GLM FIT
%         %build A matrix
%         A=zeros(N,nnz);
%         for i = 1:nnz
%             I=zeros(MM+N*ps,1);
%             I(indnz(i))=1;
%             A(:,i)=abeloperatorWS(I,1,K,B,qmf,S,divX,E,ps);
%         end
%         options.fistaitr = 1000;
%         options.fistax0=zeros(nnz,1);
%         options.indpos=1:nnz;
%         sol.glm0  = FISTA_Poi(A, O', F(:), 0, Lf , options,qmf,S,divX(indnz),0);
%         dhattemp=dhattemp*0;
%         dhattemp(indnz-MM)=sol.glm0.x./divX(indnz)';
        
        [lambdaPS, divXPS]=lambdaQUTastro_checkPS(F(:)',K,O,qmf,S,B,E);
        inddd=((MM+find(dhattemp==0)));
        divXPS(inddd)=Inf;
        A  = operator(@(x,mode) abeloperatorWS(x,mode,K,B,qmf,S,divXPS,E,ps),N,MM+N*ps);
        options=optionsF;
        options.fistax0=options.fistax0*0;
        options.fistaitr  = nitr;
        solPS  = FISTA_Poi(A, O', F(:), lambdaPS, Lf , options,qmf,S,divXPS,lambdaMax);
        sol.dhatpos=solPS.x((MM+1):end)./divXPS((MM+1):end)';
    end
        
end