function [K,B,qmf]=setUpOperators(M,N,R,bluralpha,blurR0,blurthresh,wavetype,wavepar,Kopt)
    %psf  = (1+(d/R0)^2).^(-alpha), d is distance to center of mask
    %

    %Blurring operator
    if bluralpha~=1
        R0=blurR0;
        alpha=bluralpha;
        
        %psf  = (1+(d/R0)^2).^(-alpha), d is distance to center of mask
        
        masksize=ceil(sqrt((blurthresh^(-1/alpha)-1)*R0^2));
        
        ii=[masksize:-1:0 1:masksize];
        xx=repmat(ii,length(ii),1);
        yy=repmat(ii',1,length(ii));
        
        d=sqrt(xx.^2+yy.^2);

        B=(1+(d/R0).^2).^(-alpha);
        B=B/sum(sum(B));
            
    else
        B=1;
    end
    
    %Abel operator
    if Kopt==1
        [K,x]=abelmatrix2Half2D(M,N,R,1,1);
    else
        K=0;
    end
    
    %Wavelet filter
    if isempty(wavetype) 
        wavetype='Haar';
        wavepar=1;
    end
    qmf=MakeONFilter(wavetype,wavepar);
        
end