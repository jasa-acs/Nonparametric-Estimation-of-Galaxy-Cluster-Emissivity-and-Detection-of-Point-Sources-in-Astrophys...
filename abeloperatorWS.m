function F=abeloperatorWS(x,mode,K,B,qmf,S,divX,E,ps)
    
    N=size(K,1);
    n=sqrt(N);

    M=size(K,2);
    
    %K=sparse(K);
    Wxw=0;
    Wxs=0;

    %destandardizing
  
    if mode==1
        x=x./divX';
        
        
        %BAWx
        %wavelets and splines
        if(~isempty(qmf))
            xw=x(1:M);
            x(1:M)=[];
            Wxw=IWT_PO(xw',0,qmf); %Inverse Wavelet Transform
        end
        if(~isempty(S))
            xs=x(1:M);
            x(1:M)=[];
            Wxs=(S*xs)'; %Inverse splines
        end
        Wx=Wxw+Wxs;
        
        if ps 
            d=x;
        end
        
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
        
        F=[];
        
        %Wavelet Transpose
        if(~isempty(qmf))
            WtAtBtxw=FWT_PO(AtBtx',0,qmf);
            F=[F;WtAtBtxw'];
        end

        %Spline Transpose
        if(~isempty(S))
            WtAtBtxs=(S'*AtBtx)';
            F=[F;WtAtBtxs'];
        end

        %Point sources
        if ps
            d=Btx;
            F=[F;d];
        end

        %standardize
        F=F./divX';
    end

end
    
    


