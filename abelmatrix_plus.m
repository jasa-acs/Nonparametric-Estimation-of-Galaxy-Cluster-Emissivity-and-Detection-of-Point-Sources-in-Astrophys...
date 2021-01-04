function [K, x]=abelmatrix_plus(M,R,kfactor,z,N,r,rf,firstPoint)
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