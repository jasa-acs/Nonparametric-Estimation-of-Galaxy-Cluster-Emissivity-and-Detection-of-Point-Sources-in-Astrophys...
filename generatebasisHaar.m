function S=generatebasisHaar(M)
    qmf=MakeONFilter('Haar',1);
    S=zeros(M,M);
    for i=0:((log(M)/log(2)-1))
        for j=0:((2^i)-1)
            x=zeros(1,M);x(j+1)=1;
            f=IWT_PO(x,i,qmf);
            S(2^i+j,:)=f;
        end
    end
    S(M,:)=zeros(1,M);
    S(M,(M/2):(M/2+1))=1;
    S(M,:)=S(M,:)*sum(S(M-1,:)/2);
    %S=[ones(1,M);S(1:(M-1),:)];
    %S([1 M/2],:)=S([M/2 1],:);
    S=eye(M);

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