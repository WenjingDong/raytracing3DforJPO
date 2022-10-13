clear;
close all;
U0 = 0.0468;
format long;
set(0,'DefaultLineLineWidth',2);
files = dir('../data/charney/std_*.mat');
nfile = length(files);
Ro = [0.025, 0.05, 0.1];

figure('unit','centimeter','position',[12,12,24,8]);
Uvec = U0*[1/4, 1/2,1];
subplot(1,2,1);
files = dir('../data/charney/cg_Ro*.mat');
for ii=1:length(files);
	load(['../data/charney/' files(ii).name]);
	loglog(t*Ro(ii)^2, Uvec(ii)./cg_rms);
	hold on;
end
ylim([0.05,10])
plot([0.1,0.1],ylim,'k--');
set(gca,'fontsize',12);
xlabel('$\mathrm{Ro}^2ft$','interpreter','latex','fontsize',14);
ylabel('$\overline{\varepsilon}$','interpreter','latex','fontsize',14);
lg = legend('$\mathrm{Ro} = 0.0025$', '$\mathrm{Ro} = 0.05$', '$\mathrm{Ro} = 0.1$','$\mathrm{Ro}^2ft=0.1$');
set(lg, 'interpreter', 'latex', 'fontsize', 10, 'location', 'northwest');
grid on;
xlim([1e-2,10]);
text(0.2e-2,10,'(a)','fontsize',14);
set(gca,'fontsize',12);

subplot(1,2,2);
files = dir('../data/charney/hist_N100f2Ro*.mat');
for ii=1:3
	load(['../data/charney/' files(ii).name]);
        time = round(0.1./Ro(ii).^2);
        [~,ind] = min(abs(time-t));
	loglog(BinLoc, BinLoc.*mean(Count(:,ind-10:ind+10),2));
	hold on;
end
grid on;
xlim([1,100]);
xlabel('$\omega/f$','interpreter','latex','fontsize',14);
ylabel('$e(\omega)$','interpreter','latex','fontsize',14);
tl = title('$\mathrm{Ro}^2ft = 0.1$');
set(tl,'interpreter','latex','fontsize',14);
text(0.4,10,'(b)','fontsize',14);
ylim([1e-4, 10]);
set(gca,'fontsize',12);

print -f1 -depsc -r300 ../../source/figures/charney_freqspread


