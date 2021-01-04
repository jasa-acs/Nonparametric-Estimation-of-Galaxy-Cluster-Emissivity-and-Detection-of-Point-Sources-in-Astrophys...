function [K, x]=abelmatrix2Half2D(M,N,R,kfactor,firstPoint,diagonal)
	%N is the number of pixels of F, nxn 
    %M is a power of 2, number of wavelets
	
	[K2, x]=abelmatrixFull2D(M/2,N,R,kfactor,firstPoint,diagonal);
	K1=K2(1:(size(K2,1)/2),size(K2,2):-1:1);
    K2(1:(size(K2,1)/2),:)=[];
	K2=K2;
	K=[[K1 K1*0];[K2*0 K2]];

end

%abel transformation using diagonal radius
function [K, x]=abelmatrix_diag(M,R,kfactor,z,N,r,rf,firstPoint)
	%M -- size of matrix

	N=floor(N*sqrt(2));
    R=R*sqrt(2);
    
	x=linspace(firstPoint,M-1+firstPoint,M);
	if isempty(rf)
        if kfactor==1 
            x=x*R/M;
        else
            x=(kfactor^x-1)/(kfactor^M-1)*R;
        end
		delta=R/M;
	else
		x=rf;
	end
    
    K=zeros(N,M);

    for i= 1:N
        if isempty(r) 
            y=(i-1)*R/N;
        else
            y=r(i);
        end
		s=sqrt(y^2+z^2);
        
        for j = find(x>s)

            if j==1 
                ri0=s;
            elseif x(j-1)<s
                ri0=s;
            else
                ri0=x(j-1);
            end
            
            ri1=x(j);
            
            if x(j)<s
                ri1=s;
            end
			
			K(i,j)=2*((ri1^2-s^2)^0.5-(ri0^2-s^2)^0.5);

        end
    end
    K=K(1:M,:);
end

%Abel matrix 1D 
function [K, x]=abelmatrix(M,R,kfactor,z,N,r,rf,firstPoint)
	%M -- size of matrix
	%R -- length of the radius
	if isempty(N)
        N=M;
	end

	x=linspace(firstPoint,M-1+firstPoint,M);
	if isempty(rf)
        if kfactor==1 
            x=x*R/M;
        else
            x=(kfactor^x-1)/(kfactor^M-1)*R;
        end
		delta=R/M;
	else
		x=rf;
	end
    
    K=zeros(N,M);

    for i= 1:N
        if isempty(r) 
            y=(i-1)*R/N;
        else
            y=r(i);
        end
		s=sqrt(y^2+z^2);
        
        for j = find(x>s)

            if j==1 
                ri0=s;
            elseif x(j-1)<s
                ri0=s;
            else
                ri0=x(j-1);
            end
            
            ri1=x(j);
            
            if x(j)<s
                ri1=s;
            end
			
			K(i,j)=2*((ri1^2-s^2)^0.5-(ri0^2-s^2)^0.5);

        end
    end
end


%FUNCTION: Abel Matrix for 2D 
function [K, x]=abelmatrix2D(M,R,kfactor,N,firstPoint,diagonal)
    x=[];
    K=[];
    if isempty(N) 
        N=M^2;
    end
    
    n=sqrt(N);
    
    if      diagonal==1
        Kdataname=strcat('Kdiag_M',num2str(M),'_N',num2str(N),'_R',num2str(R),'_kfactor',num2str(kfactor),'fp',num2str(firstPoint),'.mat');
    elseif  diagonal==0
        Kdataname=strcat('K_M',num2str(M),'_N',num2str(N),'_R',num2str(R),'_kfactor',num2str(kfactor),'fp',num2str(firstPoint),'.mat');
    end
                
	if exist(Kdataname,'file')
		load(Kdataname);
	else
		K=[];
		
        for i = 1:n
			z=(i-1)*R/n;
            if      diagonal==1
                [curK, x]=abelmatrix_diag(M,R,kfactor,z,n,[],[],firstPoint);
            elseif  diagonal==0
                [curK, x]=abelmatrix(M,R,kfactor,z,n,[],[],firstPoint);
            end
			K=[curK;K];
        end
		save(Kdataname, '-v7.3','K','x');
	end
end

%FUNCTION: Abel matrix 2D between (-R,-R)  and (R,R)
function [K, x] = abelmatrixFull2D(M,N,R,kfactor,firstPoint,diagonal)
	%N is the length of F,(number of pixels in the image)

	n=sqrt(N);
    
    if diagonal==1
        Kdataname=strcat('KFULLdiag_M',num2str(M),'_N',num2str(N),'_R',num2str(R),'_kfactor',num2str(kfactor),'fp',num2str(firstPoint),'.mat');
    elseif diagonal==0
        Kdataname=strcat('KFULL_M',num2str(M),'_N',num2str(N),'_R',num2str(R),'_kfactor',num2str(kfactor),'fp',num2str(firstPoint),'.mat');
    end
    
	if exist(Kdataname,'file')
		load(Kdataname);
	else
		[K2D, x]=abelmatrix2D(M,R,kfactor,(n/2)^2,firstPoint,diagonal);
		
		steps=1:(n/2):((n/2)^2);
        
		Ktop=[];
		Kbottom=[];
        for i = steps

			K2=K2D(i:(i+n/2-1),:);
			K1=K2(size(K2,1):-1:1,:);

			Ktop=[Ktop;K1;K2];
			Kbottom=[K1;K2;Kbottom];
			
        end
		K=[Ktop;Kbottom];
		
		save(Kdataname,  '-v7.3','K','x');

	end
end