
function S=generatebasis1(M,sshrink,intercept)
    %mm=10.^((-sshrink:(M-sshrink+1)));
    %mm=10.^linspace((-sshrink),100,M);
    S=zeros(2*M,M);    
    r=1:(M/2);

    b0=5;
    r0=5;
    bi=linspace(1/6,10,sqrt(2*M));
    ri=linspace(1,M/2,sqrt(2*M));
    count=0;
    for i=1:sqrt(2*M)
        for j=1:sqrt(2*M)
            count=count+1;
            f=(1+(r/ri(j)).^2).^(-3*bi(i));
            
            f=[fliplr(f) f];
            f=f-min(f);
            f=f/sum(f);
            %f=f/max(f);
            %f=f-mean(f);
            %f=f/sqrt(sum(f.^2));
            
            S(count,:)=f;
        end
    end
    
    %Assymetric King transformation
    S((M/2+1):end,1:(M/2))=0;
    S(1:(M/2),(M/2+1):end)=0;
    
    
    if(intercept)
        S=[ones(1,M);S(1:(M-1),:)];
    end
    notnanindex=find(sum(isnan(S'))==0&sum(S')~=0);
    sindex=floor(linspace(1,length(notnanindex),M));
    S=S(notnanindex(sindex),:);
    
    
    %S([1 M/2],:)=S([M/2 1],:);
    %hold off
end