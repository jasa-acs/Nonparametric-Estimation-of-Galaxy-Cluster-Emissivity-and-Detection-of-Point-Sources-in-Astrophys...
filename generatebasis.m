function S=generatebasis(M,sshrink,intercept)
    %mm=10.^((-sshrink:(M-sshrink+1)));
    %mm=10.^linspace((-sshrink),100,M);
    S=zeros(M,M);    
    r=1:(M/2);

    ri=M/2;
    b0=log(M)+0.5;  
    bi=10.^linspace(b0,-5,M);
    
    for i=1:M
        %sf=mm(i);
        %f=generatemother(M,sf);
        f=(1+(r/ri).^2).^(-3*bi(i));
        %loglog(f)
        %hold on
        f=[fliplr(f) f];
        f=f-mean(f);
        f=f/sqrt(sum(f^2));
        
        S(i,:)=f;

    end
    
    if(intercept)
        S=[ones(1,M);S(1:(M-1),:)];
    end
    %S([1 M/2],:)=S([M/2 1],:);
    %hold off
end



function f=generatemother(M,sf)

    load truef.mat
    cutf=900;
    truef=truef(1:cutf);
    x=linspace(0,1,(M/2));
    if sf ==1 
        xf=linspace(0,1,cutf);
    else
        xf=log(linspace(1,sf,cutf))/log(sf);
    end
    f2=interp1(xf,truef,x,'linear','extrap');
    mu0=[fliplr(f2) f2];
    f=mu0/sum(mu0);
    %f=f-mean(f);
    
end