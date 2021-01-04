
%Get tha value of Lambda QUT
function [lambda, divX, lambdaMax, divX1]=lambdaQUTastroPS(F,K,prevfw,O,qmf,B,E,ps)

    MC=100;
    N=length(F);
    M=size(K,2);

    %quantiles
    alphalevelPF=1/(sqrt(pi*log(M)));
    alphalevelPS=1/(M^2);
    alphalevel=1/(sqrt(pi*log(M)));
    
    indfw=find(prevfw~=0);
    u=zeros(N,length(indfw));
    I=zeros(M,1);I(1)=1;
    for i=1:indfw
        u(:,i)=abeloperatorPS(I(:,indfw(i)),1,K,B,qmf,ones(1,M),E,0);
        clear I
    end
    muhat=O'+u*prevfw(indfw);
    
        
    
    Muhat=repmat(muhat,1,MC);
    
    y=poissrnd(Muhat);
    
    %     
    %       for i=1:MC
    %          yMC=y(:,i);
    %          beta0MC=glmfit(u,yMC,'poisson','constant','off','link','identity','offset',O);
    %          Muhat(:,i)=O+u*beta0MC;
    %          %Muhat((Muhat(:,i)-O)==0,i)=1;
    %       end

    Y_muhat=(y-Muhat)./Muhat;
    clear Muhat
    clear y
    divX=ones(1,M);

      
%     %     %First step
%             bp=zeros(M,MC);
%             for i=1:MC
%                 bp(:,i)=abs(abeloperatorPS(Y_muhat(:,i),2,K,B,qmf,ones(1,M),E,0));
%             end
%     
%             divX1=quantile(bp',1-alphalevelPF);
%             %divX=divX1;
%             
%             clear bp
%             %with point sources
%     
%             if ps
%                 bpPS=zeros(N,MC);
%     
%                 for i=1:MC
%                     curbpPS=abs(abeloperatorPS(Y_muhat(:,i),2,K*0,B,qmf,ones(1,N+M),E,1));
%                     bpPS(:,i)=curbpPS((M+1):(M+N));
%                 end
%                 clear curbpPS
%                 divI1=quantile(bpPS',1-alphalevelPF);
%                 clear bpPS
%                 imind=reshape(1:N,sqrt(N),sqrt(N));
%                 nops=imind((floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)),(floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)));
%                 
%                 clear imind
%                 nops=nops(:);
%                 divI1(nops)=Inf;
%                 clear nops;
%                 divI1(E==0)=Inf;
%                 %divI=divI1;
%     
%             end
%             
    %Second step
    bp=zeros(M,MC);

    for i=1:MC     
        bp(:,i)=abs(abeloperatorPS(Y_muhat(:,i),2,K,B,qmf,[Inf divX(2:size(divX,2))],E,0));
    end

    resultsMC=max(bp);
    clear bp

    %             %%%%%
    %             %figure
    % %             X1=K*0;
    % %             labcell=cell(1,M)
    % %             for i =1:M
    % %                 X1(:,i)=abeloperatorPS(I(:,i),1,K,B,qmf,ones(1,M),E,0);
    % %                 labcell{i}='';
    % %             end
    % %             divSd=std(X1);
    % %             h=subplot(1,3,1); boxplot(bp1','labels',labcell);
    % %             title('without rescaling')
    % %             set(h,'XtickLabel',[],'YtickLabel',[]);
    % %             
    % %             subplot(1,3,3); boxplot(bp','labels',labcell);
    % %             title('quantile rescaling')
    % %             for i=1:MC     
    % %                 bp(:,i)=abs(abeloperatorPS(Y_muhat(:,i),2,K,B,qmf,[Inf divSd(2:size(divSd,2))],E,0));
    % %             end
    % %             subplot(1,3,2); boxplot(bp','labels',labcell);
    % %             title('rho-standardization')
    % %             asdsad
    %             %%%%%
            

    if ps
        divI=ones(1,N);
        imind=reshape(1:N,sqrt(N),sqrt(N));
        nops=imind((floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)),(floor(sqrt(N)/2)-floor(sqrt(N)/8)+1):(floor(sqrt(N)/2)+floor(sqrt(N)/8)));
        nops=nops(:);
        divI(nops)=Inf;
        divI(E==0)=Inf;
                 
        bpPS=zeros(N,MC);
         
        for i=1:MC
            curbpPS=abs(abeloperatorPS(Y_muhat(:,i),2,K*0,B,qmf,[ones(1,M) divI],E,1));
            bpPS(:,i)=curbpPS((M+1):(M+N));
        end

        resultsMCPS=[max(bpPS)];
        clear bpPS

        lambda=quantile(resultsMC, 1-alphalevelPF);
        lambdaI=quantile(resultsMCPS, 1-alphalevelPS);

        pd = fitdist(resultsMCPS','ExtremeValue');
        lambdaIfit = icdf(pd,1-alphalevelPS);

        [lambda lambdaI lambdaIfit]
        lambdaI=lambdaIfit;

        %standardize with pointsources
        divX=[divX divI*lambdaI/lambda];
        divX1=divX;
        divX1=[divX1 divI1*lambdaI/lambda];

    else
        lambda=quantile(resultsMC, 1-alphalevel);
        divX1=divX;
    end

    
    lambdaMax=max(abs(abeloperatorPS((F'-muhat)./muhat,2,K,B,qmf,[Inf divX(2:size(divX,2))],E,ps)));
    lambda
    %lambda=lambdaMax*0.90

end