function plotwithci(x,fmedian,fmin,fmax,fhat,f,fdom,fdommedian,fdommin,fdommax)
    M=length(fmedian);
%axis([XMIN XMAX YMIN YMAX])
    fill([x,fliplr(x)],[(fdommin),fliplr(fdommax)],[0 0 1],'EdgeColor','None','FaceAlpha',0.3)
    hold on
    fill([x,fliplr(x)],[(fmin),fliplr(fmax)],[1 0 0],'EdgeColor','None','FaceAlpha',0.3)
    axis tight
    semilogy(x,f,'color','black')
    plot(x,fdom,'color','blue','LineWidth',2)
    plot(x,fhat,'color','red','LineWidth',2)
    set(gca,'YScale','log')
    hold off
    
end