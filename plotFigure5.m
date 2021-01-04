%Code to plot the real data and compare results

XMM.pixel = 0.000694444444444445; % degree
Chandra.pixel = 0.00027333333333334; %degree
ktelescope=1500;


load('data_and_results/resultsRealData.mat')
load('data_and_results/realdata.mat')

%%%%%CHANDRA%%%%%%%%%%%5
pixel=Chandra.pixel;
Freal=realdata(2).F;
Freal=Freal(9:520,9:520);

bigF=ftel(2).Fhatfull;
bigFdom=ftel(2).Fdomfull;
fdom=ftel(2).fdomfull;

xf=linspace(-256,256,512);
x=xf*pixel;

fmedian=median(bigF);
fmax=quantile(bigF,0.975);
fmin=quantile(bigF,0.025);
fdommedian=median(bigFdom);
fdommax=quantile(bigFdom,0.975);
fdommin=quantile(bigFdom,0.025);

x=xf*pixel;
x1=x;
subplot(3,2,1)
imagesc(x,x,log(Freal))

title('Chandra')
ylabel('X-ray emission')
axis([x1(1) x1(end) x1(1) x1(end)])

subplot(3,2,3)

fill([x1,fliplr(x1)],[(ktelescope*fdommin),fliplr(ktelescope*fdommax)],[0 0 1],'EdgeColor','None','FaceAlpha',0.3)
hold on
fill([x1,fliplr(x1)],[(ktelescope*fmin),fliplr(ktelescope*fmax)],[1 0 0],'EdgeColor','None','FaceAlpha',0.3)
plot (x1,ktelescope*fdommedian,'color','blue','LineWidth',2)
plot(x1,ktelescope*fmedian,'color','red','LineWidth',2)

ylabel('Estimated log-profile')
xlabel('r')
set(gca,'YScale','log')
hold off
axis tight

fchandra=fmedian;
fchandradom=fdom;
xchandra=x1;

%%%%%%XMM%%%%%%%%%%%%%%
pixel=XMM.pixel; 
Freal=realdata(1).F;
Freal=Freal(9:520,9:520);

bigF=ftel(1).Fhatfull;
bigFdom=ftel(1).Fdomfull;
fdom=ftel(1).fdomfull;
x=xf*pixel;

fmedian=median(bigF);
fmax=quantile(bigF,0.975);
fmin=quantile(bigF,0.025);
fdommedian=median(bigFdom);
fdommax=quantile(bigFdom,0.975);
fdommin=quantile(bigFdom,0.025);
xf1=xf(abs(x)<0.16);

x1=x;
subplot(3,2,2)
imagesc(x,x,log(Freal))
title('XMM')
axis([x1(1) x1(end) x1(1) x1(end)])

subplot(3,2,4)
fill([x1,fliplr(x1)],[(fdommin),fliplr(fdommax)],[0 0 1],'EdgeColor','None','FaceAlpha',0.3)
hold on
fill([x1,fliplr(x1)],[(fmin),fliplr(fmax)],[1 0 0],'EdgeColor','None','FaceAlpha',0.3)
plot (x1,fdommedian,'--','color','blue','LineWidth',2)
plot(x1,fmedian,'--','color','red','LineWidth',2)

xlabel('r')
set(gca,'YScale','log')
hold off
axis tight

fxmm=fmedian;
fxmmdom=fdom;
xxmm=x1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%subplot(4,2,[5 6 7 8])
subplot(3,2,[5 6])
plot (xchandra,ktelescope*fchandra,'color','red','LineWidth',2)
hold on
plot (xchandra,ktelescope*fchandradom,'color','blue','LineWidth',2)

plot (xxmm,fxmm,'--','color','red','LineWidth',2)
plot (xxmm,fxmmdom,'--','color','blue','LineWidth',2)

title ('Chandra & XMM')
ylabel('Estimated log-profile')
xlabel('r')
set(gca,'YScale','log')
hold off
axis tight

legend('ALIAS - Chandra','Conventional - Chandra','ALIAS - XMM','Conventional - XMM')












