
%Get tha value of Lambda QUT
function [lambda, divX]=lambdaQUTastro_checkPS(F,K,O,qmf,S,B,E)

    MC=1000;
    N=length(F);
    M=size(K,2);
    MM=M*(~isempty(qmf))+M*(~isempty(S));
    
    %quantiles
    alphalevelPS=1/N;
    
    muhat=O';
    
    clear O
    Muhat=repmat(muhat,1,MC);
    y=poissrnd(Muhat);
    Y_muhat=(y-Muhat)./Muhat;

    clear Muhat
    clear y
    clear muhat
    
    %First step
    bpPS=zeros(N,MC);

    for i=1:MC
        curbpPS=abs(abeloperatorWS(Y_muhat(:,i),2,K*0,B,qmf,S,ones(1,N+MM),E,1));
        bpPS(:,i)=curbpPS((MM+1):(MM+N));
    end
    clear curbpPS
    divI1=quantile(bpPS',1-alphalevelPS);
    clear bpPS
    
    imind=reshape(1:N,sqrt(N),sqrt(N));
    nops=imind((floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)),(floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)));
    clear imind
    nops=nops(:);
    divI1(nops)=Inf;
    clear nops;
    
    divI1(E==0)=Inf;
    divI1(F<=1)=Inf;
    divI=divI1; 'standardizing point sources';

            
    %Second step
    bpPS=zeros(N,MC);

    for i=1:MC
        curbpPS=abs(abeloperatorWS(Y_muhat(:,i),2,K*0,B,qmf,S,[ones(1,MM) divI],E,1));
        bpPS(:,i)=curbpPS((MM+1):(MM+N));
    end

    resultsMCPS=[max(bpPS)];
    clear bpPS

    lambdaI=quantile(resultsMCPS, 1-alphalevelPS)

    
    if(sum(resultsMCPS)==MC)
        lambda=1;
    else
        pd = fitdist(resultsMCPS','ExtremeValue');
        lambdaIfit = icdf(pd,1-alphalevelPS);
        lambda=lambdaIfit
    end
    divX=[ones(1,MM)*Inf divI];
    
    if(isnan(lambda))
        warning('lambda is not numeric, fixing to the intercept')
        lambda=0;
        divX(2:end)=Inf;
        divX(1)=1;
        divX1=divX;
    end


end