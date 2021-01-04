%Code to plot the real data and compare results

XMM.pixel = 0.000694444444444445; % degree
Chandra.pixel = 0.00027333333333334; %degree


load('data_and_results/resultsRealData.mat')
load('data_and_results/realdata.mat')
ii=[2 1 3 4 5]

for i=1:5
    curi=ii(i)
    name=ftel(curi).name;
    if strcmp(ftel(curi).telescope,'chandra')
        pixel=Chandra.pixel;
        tel='Chandra';
    elseif strcmp(ftel(curi).telescope,'xmm')
        pixel=XMM.pixel;
        tel='XMM';
    end
    
    if strcmp(name,'abell2142')
        name=strcat(name,' ',tel);
    end
    
    Freal=realdata(curi).F;
    Freal=Freal(9:520,9:520);
    xf=linspace(-256,256,512);
    x=xf*pixel;
    x1=x;
    
    %Plot emission
    ax1=subplot(5,3,3*(i-1)+1)
    logFreal=log(Freal+1);
    logFreal=logFreal-min(min(logFreal));
    logFreal=floor(logFreal*300/max(max(logFreal)));
    image(logFreal)
    cm=colormap(ax1);
    logFrealcm=ind2rgb(logFreal,cm);
    imshow(logFrealcm,'XData',x,'YData',x)
    
    axis on
    
    if(i==1) 
        title('X-ray emission')
    end
    if i==5
        xlabel('r')
    end
    ylabel(name)
    axis([x1(1) x1(end) x1(1) x1(end)])
    clear ax1
    clear cm

end

if 1==2
    for i=1:5 
        curi=ii(i)
        name=ftel(curi).name;

        Freal=realdata(curi).F;
        Freal=Freal(9:520,9:520);

        Dhat=ftel(curi).Dhatfull_wavelet;

        xf=linspace(-256,256,512);
        x=xf*pixel;
        x1=x;


        %Plot Point source detection
        ax2=subplot(5,3,3*(i-1)+2)

        imshow(log(Freal+1)/max(max(log(Freal+1))));

        dhatcount=sum(Dhat>0);
        %ind=find(dhatcount>(size(Dhat,1)/5));

        dhat=max(Dhat);
        %thresh = multithresh(dhat(dhat~=0),1)
        %ind=find(dhat>thresh(end)&dhatcount>(size(Dhat,1)/5));
        ind=find(dhat>max(dhat/2));
        ind=find(dhat>quantile(dhat(dhat~=0),0.75));
        %thresh=mean(dhat(dhat~=0))+mad(dhat(dhat~=0));
        %ind=find(dhat>thresh);

        indX=ceil(ind/512);
        indY=mod(ind,512);
        hold on
        scatter(indX,indY,3)
        hold off
        if i==1 
            title('Residuals')
        end
        if i==5
            xlabel('r')
        end

    end
else
    for i=1:5 
        curi=ii(i)
        name=ftel(curi).name;
        Freal=realdata(curi).F;
        Freal=Freal(9:520,9:520);
        
            name=ftel(curi).name;
    if strcmp(ftel(curi).telescope,'chandra')
        pixel=Chandra.pixel;
        tel='Chandra';
    elseif strcmp(ftel(curi).telescope,'xmm')
        pixel=XMM.pixel;
        tel='XMM';
    end
        
        Ereal=realdata(curi).E;
        Ereal=Ereal(9:520,9:520);
        
        indmax=find(Freal==max(max(Freal)),1);

        Dhat=ftel(curi).Dhatfull_wavelet;
        Dhat2=reshape(max(Dhat),512,512);
        
        %Dhat2(indmax)=Freal(indmax);
        Dhat2(Dhat2>2*Freal(256,256))=2*Freal(256,256);
        Dhat2(256,256)=2*Freal(256,256);
        
        xf=linspace(-256,256,512);
        x=xf*pixel;
        x1=x;


        %Plot Point source detection
        ax2=subplot(5,3,3*(i-1)+2)
        dhat=max(Dhat);
        %dhat=dhat/max(dhat);
        logdhat=log(dhat+1);
        %logdhat=logdhat/max(logdhat);
        ind=find(dhat>max(dhat/2));
        ind=find(dhat>quantile(dhat(dhat~=0),0.75));

        dhatnew=dhat*0;
        dhatnew(ind)=dhat(ind);
        
        dhatnew=reshape(dhatnew,512,512);
        
        axis off
        %imagesc(log(dhatnew))
        
        logDreal=log(dhatnew+1);
        logDreal=logDreal-min(min(logDreal));
        logDreal=floor(logDreal*300/max(max(logDreal)));
        %image(logDreal,'XData',x,'YData',x)
        %image(reshape(logdhat*7,512,512))
        %plot(logdhat)
        logdhat2=logdhat*0;
        logdhat2(ind)=logdhat(ind);
        logdhat2=logdhat2-min(logdhat2(ind));
        logdhat2(logdhat2<0)=0;
        %logdhat(ind)=max(logdhat);
        %image(reshape(logdhat,512,512))

        %imshow(reshape(dhat/max(dhat),512,512))
        %mesh(reshape(dhat,512,512))
        %imagesc(Dhat2);
        
        %axis on
        %axis([x1(1) x1(end) x1(1) x1(end)])
        DDD=1-Dhat2/(max(max(Dhat2)));
        %DDD(1,:)=0;
        %DDD(512,:)=0;
        %DDD(:,1)=0;
        %DDD(:,512)=0;
        
        imshow(DDD,'XData',x,'YData',x)
        axis on
        if i==1 
            title({'$\hat{S}$'},'Interpreter','latex')
        end
        if i==5
            xlabel('r')
        end

    end
end


for i=1:5 
    curi=ii(i)
    name=ftel(curi).name;
    ktelescope=1;
    if strcmp(ftel(curi).telescope,'chandra')
        pixel=Chandra.pixel;
        tel='Chandra';
        ktelescope=1500;
    elseif strcmp(ftel(curi).telescope,'xmm')
        pixel=XMM.pixel;
        tel='XMM';
    end
    
    if strcmp(name,'abell2142')
        name=strcat(name,' ',tel);
    end
    
    xf=linspace(-256,256,512);
    x=xf*pixel;
    x1=x;
    
    bigF=ftel(curi).Fhatfull;
    bigFdom=ftel(curi).Fdomfull;
    fdom=ftel(curi).fdomfull;
    fmedian=median(bigF);
    fmax=quantile(bigF,0.975);
    fmin=quantile(bigF,0.025);
    fdommedian=median(bigFdom);
    fdommax=quantile(bigFdom,0.975);
    fdommin=quantile(bigFdom,0.025);
    
    %Plot estimated profile with confidence intervals
    ax3=subplot(5,3,3*(i-1)+3)
    fill([x1,fliplr(x1)],[(ktelescope*fdommin),fliplr(ktelescope*fdommax)],[0 0 1],'EdgeColor','None','FaceAlpha',0.3)
    hold on
    fill([x1,fliplr(x1)],[(ktelescope*fmin),fliplr(ktelescope*fmax)],[1 0 0],'EdgeColor','None','FaceAlpha',0.3)
    plot (x1,ktelescope*fdommedian,'color','blue','LineWidth',2)
    plot(x1,ktelescope*fmedian,'color','red','LineWidth',2)

    if i==1 
        title('Estimated log-profile')
    end
    if i==5
        xlabel('r')
    end
    
    
    set(gca,'YScale','log')
    hold off
    axis tight
    clear ax2
    clear ax3
end










