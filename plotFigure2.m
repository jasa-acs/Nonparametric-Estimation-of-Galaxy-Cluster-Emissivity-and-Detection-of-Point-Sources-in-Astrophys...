M=512
subplot(2,4,1)
x1=linspace(-1,1,M);
x=linspace(0,1,M/2);
load('data_and_results/rsummary_diag_powerlaw_512blockssymps1type1nboot24.mat')
semilogy(x,summary.f,'k')
axis tight
title('cosmoBlocks')
ylabel('Profile function')
xlabel('r')

subplot(2,4,5)
imagesc(log(curresult.F),'XData',x1,'YData',x1)
ylabel('log X-ray emission')
xlabel('r')

subplot(2,4,2)
load('data_and_results/rsummary_diag_powerlaw_512truefps1type1nboot24.mat')
semilogy(x,summary.f,'k')
axis tight
title('cosmo1')
xlabel('r')

subplot(2,4,6)
imagesc(log(curresult.F),'XData',x1,'YData',x1)
xlabel('r')


subplot(2,4,3)
load('data_and_results/rsummary_diag_powerlaw_512f4symps1type1nboot24.mat')
semilogy(x,summary.f,'k')
axis tight
title('cosmo2')
xlabel('r')

subplot(2,4,7)
imagesc(log(curresult.F),'XData',x1,'YData',x1)
xlabel('r')

subplot(2,4,4)
x=linspace(-1,1,M);
load('data_and_results/rsummary_diag_powerlaw_512f4ps1type1nboot24.mat')
semilogy(x,summary.ffull,'k')
axis tight
title('cosmoAsym')
xlabel('r')

subplot(2,4,8)
imagesc(log(curresult.F),'XData',x1,'YData',x1)
xlabel('r')
