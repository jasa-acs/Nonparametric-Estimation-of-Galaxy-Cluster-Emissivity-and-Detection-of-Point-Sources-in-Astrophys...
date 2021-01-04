
%Get tha value of Lambda QUT
function [lambda, divX, lambdaMax, divX1]=lambdaQUTastroWS(F,K,prevfw,O,qmf,S,B,E,ps)

    bps=[];bpw=[];
    divXw=[];divXs=[];
    MC=100;
    N=length(F);
    M=size(K,2);
    MM=M*(~isempty(qmf))+M*(~isempty(S));

    %quantiles
    alphalevelPF=1/(sqrt(pi*log(M)));
    alphalevelPS=1/(M);
    indfw=find(prevfw~=0);
    u=zeros(N,length(indfw));
    I=zeros(MM,1);I(1)=1;
    for i=1:length(indfw)
        u(:,i)=abeloperatorWS(I(:,indfw(i)),1,K,B,qmf,S,ones(1,MM),E,0);
        clear I
    end

    muhat=O'+u*prevfw(indfw);
    clear prevfw
    clear O
    clear u
    if(~isempty(qmf)) 
        indfww=indfw(indfw<=M);
        indfws=indfw(indfw>M)-M;
    else
        indfws=indfw(indfw<=M);
    end
    Muhat=repmat(muhat,1,MC);
    y=poissrnd(Muhat);
    Y_muhat=(y-Muhat)./Muhat;

    clear Muhat
    clear y
    clear muhat
    
    divX=[];
    divX1=[];
    
    %First step
        %wavelets
        if(~isempty(qmf))
            divXw=ones(1,M);
            bp=zeros(M,MC);
            for i=1:MC
                bp(:,i)=abs(abeloperatorWS(Y_muhat(:,i),2,K,B,qmf,[],ones(1,M),E,0));
            end
            
            divX1w=quantile(bp',1-alphalevelPF);
            divX1w=ones(1,M)*max(divX1w(2:end));'max rescaling wavelets'
            divXw=divX1w; 'standardizing wavelets'
            divXw((M/2+2):M)=Inf;'removing high frequency wavelets'

        end
        
        %splines
        if(~isempty(S))
            divXs=ones(1,M);
            bp=zeros(M,MC);
            for i=1:MC
                bp(:,i)=abs(abeloperatorWS(Y_muhat(:,i),2,K,B,[],S,ones(1,M),E,0));
            end

            divX1s=quantile(bp',1-alphalevelPF);
            divXs=divX1s; 'standardizing splines'
            %divXs(M)=Inf;
            divXs(2:2:M)=Inf;'removing high frequency splines'

        end

        clear bp
       
        %with point sources
        if ps
            bpPS=zeros(N,MC);

            for i=1:MC
                curbpPS=abs(abeloperatorWS(Y_muhat(:,i),2,K*0,B,qmf,S,ones(1,N+MM),E,1));
                bpPS(:,i)=curbpPS((M+1):(M+N));
            end
            clear curbpPS
            divI1=quantile(bpPS',1-alphalevelPF);
            clear bpPS
            imind=reshape(1:N,sqrt(N),sqrt(N));
            nops=imind((floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)),(floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)));

            clear imind
            nops=nops(:);
            divI1(nops)=Inf;
            clear nops;
            divI1(E==0)=Inf;
            divI=divI1; 'standardizing point sources';
        end
            
    %Second step
        bp=zeros(M,MC);
        lambdaw=0;
        %wavelets
        if(~isempty(qmf))
            for i=1:MC     
                divLambda=divXw;
                divLambda(indfww)=Inf;
                bp(:,i)=abs(abeloperatorWS(Y_muhat(:,i),2,K,B,qmf,[],divLambda,E,0));
            end
            bpw=bp;
            resultsMCw=max(bp);
            lambdaw=quantile(resultsMCw, 1-alphalevelPF)
            divX=divXw;
            divX1=divX1w;
            
            lambda=lambdaw;
        end
        
        %splines
        if(~isempty(S))
            for i=1:MC     
                divLambda=divXs;
                divLambda(indfws)=Inf;
                bp(:,i)=abs(abeloperatorWS(Y_muhat(:,i),2,K,B,[],S,divLambda,E,0));
            end
            bps=bp;
            resultsMCs=max(bp);
            lambdas=quantile(resultsMCs, 1-alphalevelPF)
            
            divX=[divX*lambdaw/lambdas divXs];
            divX1=[divX1*lambdaw/lambdas divX1s];
            
            lambda=lambdas;
        end
        
        resultsMC=max([bpw; bps]);
        lambda=quantile(resultsMC, 1-alphalevelPF)
        %lambda=quantile(resultsMC, 1-0.01),'0.01'
        divX=[divXw divXs];
        

        %clear bp

        if ps
            divI=ones(1,N);
            imind=reshape(1:N,sqrt(N),sqrt(N));
            nops=imind((floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)),(floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)));
            nops=nops(:);
            divI(nops)=Inf;
            divI(E==0)=Inf;

            bpPS=zeros(N,MC);

            for i=1:MC
                curbpPS=abs(abeloperatorWS(Y_muhat(:,i),2,K*0,B,qmf,S,[ones(1,MM) divI],E,1));
                bpPS(:,i)=curbpPS((MM+1):(MM+N));
            end

            resultsMCPS=[max(bpPS)];
            bp1=bpPS(sum(bpPS')~=0,:);
             boxplot([bp' bp1(1:1000,:)'])
             
            clear bpPS
    
            lambdaI=quantile(resultsMCPS, 1-alphalevelPS)

            pd = fitdist(resultsMCPS','ExtremeValue');
            lambdaIfit = icdf(pd,1-alphalevelPS);

            [lambdaI lambdaIfit]
            lambdaI=lambdaIfit;

            %standardize with pointsources
            divX=[divX divI*lambdaI/lambda];
            divX1=divX;
            divX1=[divX1 divI1*lambdaI/lambda];

        end

    divLambdaMax=divX;
    divLambdaMax(indfw)=Inf;
    %lambdaMax=max(abs(abeloperatorspline((F'-muhat)./muhat,2,K,B,S,[Inf divX(2:size(divX,2))],E,ps)));
    %lambdaMax=max(abs(abeloperatorWS((F'-muhat)./muhat,2,K,B,qmf,S,divLambdaMax,E,ps)));
    if(isnan(lambda))
        warning('lambda is not numeric, fixing to the intercept')
        lambda=0;
        divX(2:end)=Inf;
        divX(1)=1;
        divX1=divX;
    end
    lambdaMax=-1;
    %lambda
    %lambda=lambdaMax*0.99

end