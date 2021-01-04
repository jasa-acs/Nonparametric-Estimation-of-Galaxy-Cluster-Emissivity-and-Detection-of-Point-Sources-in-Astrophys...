function [sol,F]= astrosolvePS(F,R,E,O,K,options)
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
%K       -> K matrix (computed if length(k)==1)

%options ->Structure with the following options. Check defaults with astro_setparams()

    %bluralpha->  parameter for psf (Set to 0 for no blurring)
    %blurR0-> parameter for psf
    %blurthresh->  parameter for psf
    
    %wavelet -> 2 entries vector. First entry is the wavelet type and second is
    %           wavelet parameter.  Wavelet type is 1 for 'Daubechies', 
    %           2 for 'Symmlet', 3 for 'Coiflet'. If length(wavelet)=1 then
    %           wavelet is set to [1 8].
    %nitr    -> Number of iterations for convergence, nitr=200 is enough but
    %           nitr=1000 is better to be sure
    %O       ->	Background. It must be in vector form of size of a 1xn^2 vector. If O=0 then 
    %           there is no offset. If the background is an image Im, then you should do O=Im(:)'
    %E       ->	Sensitivity. It must be in vector form of size of a 1xn^2 vector. If E=1 then 
    %           there is no sensitivity effect. If the background is an image Im, then you should do E=Im(:)'
    %ps      -> If there are point sources in image, set to 1.  Set to 0 if not
    %psclean -> Set to 1 if point sources cleaning is desired after FISTA and GLM
    %smoothf -> If a smoothing filter should be applied after profile estimation
    %showplot-> Show plot
    %prevfw  -> value of not penalized wavelet coefficients (taken into account
    %           if length(prevfw)~=0))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
%Structure with the following elemets:
%
%logf     -> Estimated log profile
%fhat     -> Estimated profile
%fwhat    -> Estimated wavelet coefficients using GLM
%dhat     -> point sources using GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    remblack  =opts.remblack

    %Check Image size, make them squares of even size, convert to vectors
    [F,E,O]=astro_setFEO(F,E,O,1,opts.remblack);
    R=size(F,1)/2;
    %Convert images in vectors
    %F(F==0)=1e-10;
    F=F(:)';
    O=O(:)';
    E=E(:)';
    
    N=length(F);
    n=sqrt(N);
    M=2^floor(log(sqrt(N))/log(2));
    Lf=[];
    
    %weights, not used
    w=ones(1,length(F))';
    
    %set wavelet type and parameter
    if length(wavelet)<2 wavelet=[1 8]; end
    wavetype=wavelet(1);
    wavepar= wavelet(2);
    if wavelet(1)==0
       wavetype='Haar'; 
    elseif wavelet(1)==1
        wavetype='Daubechies';
    elseif wavelet(1)==2
        wavetype='Symmlet';
    elseif wavetype==3
        wavetype='Coiflet';
    end
    
    %Set up operators
    if isempty(K)
        [K,B,qmf]=setUpOperators(M,N,R,bluralpha,blurR0,blurthresh,wavetype,wavepar,1);
    else
        [K,B,qmf]=setUpOperators(M,N,R,bluralpha,blurR0,blurthresh,wavetype,wavepar,0);
    end
    
    % FISTA OPTIONS
    options.fistaitr  = nitr;
    options.fistatol  = 1e-60;
    options.fistalogperiod   = 1;
    
    %Estimate intercept
    I=zeros(M,1);    I(1)=1;
    itersolver=0;
    
    if isempty(prevfw)
        itersolver=1; %iterate solver
        A=abeloperatorPS(I,1,K,B,qmf,ones(1,M),E,0); 

        %beta0 is the estimate of the father wavelet coefficient under the null model
        if ps
            beta0=glmfit(A,F','poisson','constant','off','link','identity','offset',O,'weights',w') 
            beta0=sum(F~=0)*median(F(F~=0))/(sum(A))
        else
            beta0=glmfit(A,F','poisson','constant','off','link','identity','offset',O,'weights',w') 
        end
       
        prevfw=0*ones(M+N*ps,1);
        prevfw(1)=beta0;
    end
    
    clear A
    %Get Lambda and rescaling factors
    [lambda, divX,lambdaMax,divX1]=lambdaQUTastroPS(F,K,prevfw,O,qmf,B,E,ps);

    %SOLVE GLM L1 
    
    %Atemp=zeros(N,M);
    %for i = 1:M
    %    I=zeros(M,1);I(i)=1;
    %    Atemp(:,i)=abeloperatorPS(I,1,K,B,qmf,ones(1,M),E,0); 
    %end
    %betaini=glmfit(Atemp,F','normal','constant','off','link','identity','offset',O,'weights',w');

    %load beta0.mat
    %prevfw=beta0;
    options.fistax0=prevfw.*divX';
    options.fistax0(isnan(options.fistax0))=0;
    
    A  = operator(@(x,mode) abeloperatorPS(x,mode,K,B,qmf,divX,E,ps),N,M+N*ps);
    A1=A;
    sol.static  = FISTA_Poi(A, O', F', lambda, Lf , options);
    sol.static.x=sol.static.x./divX';

    %setup estimates
    sol.glm=sol.static.x*0; 
    fwlasso=sol.static.x(1:M);
    sol.static.x(1:M)=[];
    
    dlassopure=sol.static.x;
    dlasso2=sol.static.x;

    %Detect point sources
    if ps
        dlasso2(dlasso2<0)=0;
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
    clear dlasso2
    nnz=sum(sol.static.x~=0);
    indnz=find(sol.static.x~=0);
    
    %GET GLM FIT
        %build A matrix
        A=zeros(N,nnz);
        for i = 1:nnz
            I=zeros(M+N*ps,1);
            I(indnz(i))=1;
            A(:,i)=abeloperatorPS(I,1,K,B,qmf,divX,E,ps);
        end
        
        %get beta from glm with MATLAB
        betanew=glmfit(A,F','poisson','constant','off','link','identity','offset',O','weights',w);
        
        %if negative objective function, get with FISTA
        if imag(sum(A*betanew+O' - F' .* log(A*betanew+O')))~=0
            %starting point assuming gaussian
            betanew=glmfit(A,F','normal','constant','off','link','identity','offset',O','weights',w);
            options.fistaitr  = 1000;
            options.fistax0 = betanew;
            sol.glm0  = FISTA_Poi(A, O', F', 0, Lf , options);
            if sol.glm0.conv==0
                options.fistax0=sol.static.x(sol.static.x~=0);
                sol.glm0  = FISTA_Poi(A, O', F', 0, Lf , options);
            end
            if sol.glm0.conv==0
            	betanew=sol.static.x(sol.static.x~=0);
            else
                betanew=sol.glm0.x;   
            end
        end       
    
    %get estimates
    
    sol.glm(indnz)=betanew./divX(indnz)';
    fwhat=sol.glm(1:M);
    sol.glm(1:M)=[];
    dhat=sol.glm;
    if ps==0 dhat=zeros(1,N); end

    fhat=IWT_PO(fwhat,0,qmf);
    flasso=IWT_PO(fwlasso,0,qmf);
    clear sol;

    sol.fhat=fhat;
    sol.fwhat=fwhat;

    sol.fwlasso=fwlasso;
    sol.flasso=flasso;
    sol.dlasso=dlassopure;
    sol.dhat=dhat;

    %get log(f) and interpolate negative values of f
    y=log(fhat);
    y(fhat<0)=NaN;
    x=linspace(-M/2,M/2,M);
    xi=x(find(~isnan(y)));yi=y(find(~isnan(y)));
    logf=interp1(xi,yi,x,'linear','extrap');

    if smoothf  logf=smooth(logf); end;
    sol.logf=logf;

    sol.options=options;

    if showplot astroplot(sol,reshape(F,sqrt(N),sqrt(N))); end;
    sol.R=R;
    sol.nboot=0;
    F=reshape(F,sqrt(length(F)),sqrt(length(F)));
end