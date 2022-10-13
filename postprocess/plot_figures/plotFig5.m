clear;
close all;
set(0,'DefaultLineLineWidth',2);
figure('unit','centimeter','position',[10,10,24,8]);
subplot(1,2,1);
files = dir('../data/charney/omegaavg_init*.mat');
for ii=1:length(files)
	load(['../data/charney/' files(ii).name]);
        [tavg, omegaravg] = timeavgk2(t, omegar_avg, 10);
	plot(tavg, omegaravg);
	hold on;
end
grid on;
xlabel('$ft$','interpreter','latex','fontsize',14);
ylabel('$\bar{\omega}$','interpreter','latex','fontsize',14);
xlim([0,3000]);
ylim([1,5]);
grid on;
text(-500, 5, '(a)','fontsize', 14);
set(gca,'fontsize',12);

subplot(1,2,2);
files = dir('../data/charney/hist_init_*.mat');
for ii=1:length(files)
	load(['../data/charney/' files(ii).name]);
	loglog(BinLoc, BinLoc .* mean(Count(:,end-20:end),2));
	hold on;
end
xlabel('$\omega/f$','interpreter','latex','fontsize',14);
ylabel('$e(\omega)$','interpreter','latex','fontsize',14);
grid on;
lg = legend('$\mathrm{Ro} = 0.1, \omega_0=2f$',...
'$\mathrm{Ro} = 0.1, \omega_0=4f$',...
'$\mathrm{Ro} = 0.1, \omega_0=5f$',...
'$\mathrm{Ro} = 0.2, \omega_0=2f$');
set(lg, 'interpreter', 'latex','fontsize',12);
text(0.4,10,'(b)','fontsize',14);
set(gca,'fontsize',12);

print -f1 -depsc -r300 ../../../source/figures/charney_universal
