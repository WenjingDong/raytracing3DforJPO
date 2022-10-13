clear;
close all;
set(0,'defaultLineLinewidth',2);
figure('unit','centimeter','position',[10,10,12,8]);
files = dir('../data/QG/eigen_*');
slopes = [];
for ii=1:length(files)
	load(['../data/QG/' files(ii).name]);
	y = BinLoc.*pdfValues;
	loglog(BinLoc, BinLoc.*pdfValues);
        slope = compute_slope(BinLoc, y);
        slopes = [slopes;slope];	
	hold on;
end

files = dir('../data/QG/hist*');
for ii=1:length(files)
	load(['../data/QG' files(ii).name]);
	y = BinLoc.*mean(Count(:,1003-10:1003+10),2);
	loglog(BinLoc, y);
        slope = compute_slope(BinLoc, y); 
        slopes = [slopes;slope];	
	hold on;
end

lg = legend('QG flow, eigenvector','synthetic flow, eigenvector',...
'QG flow, ray tracing', 'synthetic flow, ray tracing');
set(gca,'fontsize',12);
xlabel('$\omega/f$','interpreter','latex','fontsize',14);
ylabel('$e(\omega)$','interpreter','latex','fontsize',14);
grid on;
print -f1 -depsc -r300 ../../../source/figures/QG_syn_stationary
