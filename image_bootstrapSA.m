function [Fdom, Fdomfull]=image_bootstrapSA(F1,E1,O1,M,center,bootwindow,nboot)
    Fdom=zeros(nboot,M/2);
    Fdomfull=zeros(nboot,M);
    for i=1:nboot
        strcat('bootstrap ',num2str(i))
        [F,E,O,~,~]=image_bootstrap(F1,E1,O1,center,bootwindow);
        M=size(F,1);
        
        dominique=astroStateOfTheArt(F,E,O,[0 360]);
        fdom=interpdom(dominique,M)';
        
        dominique=astroStateOfTheArt(F,E,O,[0 180]);
        fdomd=interpdom(dominique,M)';
        
        dominique=astroStateOfTheArt(F,E,O,[180 360]);
        fdomi=interpdom(dominique,M)';
        
        fdomfull=[fliplr(fdomi) fdomd];
        fdomfull=interp1([(-M/2):(-1) 1:(M/2)],fdomfull,linspace(-M/2,M/2,M));
        
        Fdom(i,:)=fdom;
        Fdomfull(i,:)=fdomfull;
    end
end

function fdom=interpdom(dominique,M)
    fdom=interp1(dominique.xbin,dominique.fbin,1:(M/2),'linear','extrap')';
    if sum(abs(fdom)==0)
        fdom=fdom+1e-100;
    else
        fdom(fdom<=0)=min(dominique.fbin(dominique.fbin>0));
    end
end

function [F,E,O,up_move,right_move]=image_bootstrap(F,E,O,center,bootwindow)
    n=size(F,1);
    if center>0
        newn=n-ceil(center*n);
        up_move=floor(rand(1)*ceil(center*n));
        right_move=floor(rand(1)*ceil(center*n));
        F=F((1+up_move):(up_move+newn),(1+right_move):(right_move+newn));
        E=E((1+up_move):(up_move+newn),(1+right_move):(right_move+newn));
        O=O((1+up_move):(up_move+newn),(1+right_move):(right_move+newn));
    else
        up_move=0;right_move=0;
    end
    n=size(F,1);
 
%% 2 points symmetric bootstrap
    for i=1:(n/2)
        for j=1:n
            symInd=[i,(n-i+1)];
            bootInd=randsample(symInd,2,true);
            F(symInd,j)=F(bootInd,j);
            E(symInd,j)=E(bootInd,j);
            O(symInd,j)=O(bootInd,j);
        end
    end
    

    
%% blurring mask window

%     F1=F;E1=E;O1=O; 
%     masksize=(size(bootwindow,1)-1)/2;
% 
%     for i=(masksize+1):(n-masksize)
%         for j=(masksize+1):(n-masksize)
%             i0=i-masksize;
%             i1=i+masksize;
%             j0=j-masksize;
%             j1=j+masksize;
%             
%             Fwindow=F1(i0:i1,j0:j1);
%             Ewindow=E1(i0:i1,j0:j1);
%             Owindow=O1(i0:i1,j0:j1);
% 
%             ijwindow=randsample(1:size(bootwindow,1)^2,1,true,bootwindow(:));
%             
%             F(i,j)=Fwindow(ijwindow);
%             E(i,j)=Ewindow(ijwindow);
%             O(i,j)=Owindow(ijwindow);
% 
%         end
%     end

%% nxn window bootstrap

    bootwindow=4;
    for i=1:(n/bootwindow)
        for j=1:(n/bootwindow)
            i0=bootwindow*(i-1)+1;
            i1=bootwindow*i;
            j0=bootwindow*(j-1)+1;
            j1=bootwindow*j;
            
            Fwindow=F(i0:i1,j0:j1);
            Ewindow=E(i0:i1,j0:j1);
            Owindow=O(i0:i1,j0:j1);
            ijwindow=randsample(1:bootwindow^2,bootwindow^2,true);
            
            Fwindow(:)=Fwindow(ijwindow);
            Ewindow(:)=Ewindow(ijwindow);
            Owindow(:)=Owindow(ijwindow);
            
            F(i0:i1,j0:j1)=Fwindow;
            E(i0:i1,j0:j1)=Ewindow;
            O(i0:i1,j0:j1)=Owindow;
        end
    end
    
end