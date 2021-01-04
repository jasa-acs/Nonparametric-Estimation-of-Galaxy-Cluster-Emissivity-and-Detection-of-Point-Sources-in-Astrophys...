function F=abeloperatorPS(x,mode,K,B,qmf,divX,E,ps)
    
    N=size(K,1);
    n=sqrt(N);
    M=size(K,2);

    K=sparse(K);
    
    %destandardizing
  
    if mode==1
        x=x./divX';
        if ps 
            d=x((M+1):end);
            x=x(1:M);
        end
        
        %BAWx
        M=size(x,1);
        
        %Inverse Wavelet Transform
        Wx=IWT_PO(x',0,qmf);

        %Abel Transform
        AWx=K*Wx';
        
        
        %Blurring
        if length(B)>1 
            BAWx=reshape(AWx',n,n);
            BAWx= filter2(B,BAWx);
            BAWx=BAWx(:);
        else
            BAWx=AWx;
        end
       
        F=E'.*BAWx;
        
        
        %Point sources
        if ps
            
            if length(B)>1
                d=reshape(d',n,n);
                d=filter2(B,d); 
                d=d(:);
            end
            
            d=E'.*d;
            F=F+d;
        end
  
  
    elseif mode==2
        %W'A'B'x
        M=sqrt(size(x,1));
        
        %Sensitivity
        x=E'.*x;
        
        
        %Blurring Transpose
        if length(B)>1 
            y=reshape(x',M,M);
            Btx= filter2(B,y); 
            Btx=Btx(:);
        else
            Btx=x;
        end
        
        %Abel Transpose
        AtBtx=K'*Btx;
        
        %Wavelet Transpose
        WtAtBtx=FWT_PO(AtBtx',0,qmf);
        
        F=WtAtBtx';

        %Point sources
        if ps
            d=Btx;
            F=[F;d];
        end

        %standardize
        F=F./divX';
    end

end
    
    


