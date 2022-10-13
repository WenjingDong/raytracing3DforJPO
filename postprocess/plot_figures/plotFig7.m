clear;
close;
set(0, 'DefaultLineLineWidth',2);
files = dir('../data/charney/hist_a*.mat');
figure('unit','centimeter','position', [10,5, 12, 8]);
for ii=1:length(files)
	load(['../data/charney/' files(ii).name]);
	loglog(BinLoc, BinLoc .* mean(Count(:,end-50:end),2));
	hold on;
end
load('../data/charney/eigen_hist.mat');
loglog(BinLoc, BinLoc .* pdfValues);
lg = legend('$\sigma/f = 0$', '$\sigma/f = 0.1$', '$\sigma/f = 0.2$','eigenvector analysis');
set(lg, 'interpreter', 'latex');
xlabel('$\omega/f$', 'interpreter', 'latex', 'fontsize', 14);
ylabel('$e(\omega)$', 'interpreter', 'latex', 'fontsize', 14);
set(gca, 'fontsize', 12);
grid on;
ylim([1e-3, 10]);
print -f1 -depsc -r300 ../../../source/figures/charneysigma
