function [F, E, O] = astro_setFEO(F,E,O,makepair,removeblack)

    if length(E)==1 
        E=E*ones(size(F,1),size(F,1));
    else
        if size(F,1)~=size(E,1)||size(F,2)~=size(E,2)
            error('Image and Sensitivity must be the same size')
        end
    end
    
    if length(O)==1 
        O=O*ones(size(F,1),size(F,1));
    else
        if size(F,1)~=size(O,1)||size(F,2)~=size(O,2)
            error('Image and Background must be the same size')
        end
    end
    
    if removeblack
        %Remove black borders
        Frhalf=F(ceil(size(F,1)/2),:);
        Fchalf=F(:,ceil(size(F,2)/2));
        zz=find(Frhalf~=0);
        y1=zz(1);
        y2=zz(end);
        yy=min(y1,size(F,1)-y2+1);
        zz=find(Fchalf~=0);
        x1=zz(1);
        x2=zz(end);
        xx=min(x1,size(F,1)-x2+1);

        sizeF=size(F);

        F=F(yy:(sizeF(1)-yy+1),:);
        O=O(yy:(sizeF(1)-yy+1),:);
        E=E(yy:(sizeF(1)-yy+1),:);

        E=E(:,xx:(sizeF(2)-xx+1));
        F=F(:,xx:(sizeF(2)-xx+1));
        O=O(:,xx:(sizeF(2)-xx+1));
    end
    
    %Check Image sizes
    F=make2square(F,makepair);

    %sensitivity
    E=make2square(E,makepair);
    
    %ii=[(size(F,1)/2):-1:1 1:(size(F,1)/2)];
    %xx=repmat(ii,length(ii),1);
    %yy=repmat(ii',1,length(ii));
    %d=sqrt(xx.^2+yy.^2);
    %E=E.*(d<=(size(F,1)/2));
    
    %offset
    O=make2square(O,makepair);
    O=O-min(min(O))+1e-4;
    
    
    
end

function X=make2square(X,makepair)

%Check Image size
    
    if size(X,1)~=size(X,2)
        %Crop image when it is not a square
        warning('Image is being cropped to a square')
        n1=size(X,1);
        n2=size(X,2);
        n0=min(n1,n2);
        X=X(((n1-n0)/2+1):((n1-n0)/2+n0),:);
        X=X(:,((n2-n0)/2+1):((n2-n0)/2+n0));
    end
    
    if mod(size(X,1),2)~=0&& makepair
        %Crop image when side is even
        warning('Removing last column and row of pixels to get even sides')
        X(:,end)=[];
        X(end,:)=[];
    end
 
end
    
    
    