clear;
close all;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesTitleFontWeight','normal');
files = dir('velocity_*.mat');
nfile = length(files);
N = 100;
f = 1;
legendlabels = {'QG', 'random phase'};
figure('unit','centimeter','position',[10,10,14,8]);
for ii=1:nfile
	load([files(ii).name]);
	[BinLoc, pdfValues, omegarFilt, DiscMatFilt] = eigen_hist(f, N, ux, uy, uz, vx, -ux, vz);
	loglog(BinLoc, BinLoc.*pdfValues);
	hold on;
end
legend(legendlabels);
grid on;
ylabel('energy density $e(\omega_r)$','interpreter','latex','fontsize',12);
%ylim([0.1,1]);
xlabel('$\omega_r/f$','interpreter','latex','fontsize',12);
return;













