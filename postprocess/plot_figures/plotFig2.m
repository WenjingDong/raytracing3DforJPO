clear;
close all;
set(0,'defaultlinelinewidth',2);
set(0,'defaultaxestitlefontweight','normal');

figure('unit','centimeter','position',[10,10, 27,8]);
load('../data/QG/hist_QG.mat');
times = [101, 251, 602, 802, 1003];
subplot(1,2,1);
loglog(BinLoc, BinLoc .* mean(Count(:,1),2),'*');
hold on;
for ii=1:length(times)
    loglog(BinLoc, BinLoc .* mean(Count(:, times(ii)-10:times(ii)+10),2));
    hold on;
end
labels = {'$ft$ = 0 ',...
['$ft$ = ' num2str(round(t(times(1))))],...
['$ft$ = ' num2str(round(t(times(2))))],...
['$ft$ = ' num2str(round(t(times(3))))],...
['$ft$ = ' num2str(round(t(times(4))))],...
['$ft$ = ' num2str(round(t(times(5))))]};
lg = legend(labels);
set(lg,'interpreter','latex');
grid on;
xlabel('$\omega/f$','interpreter','latex','fontsize',14);
ylabel('$e(\omega)$','interpreter','latex','fontsize',14);
ylim([1e-4, 1e2]);
slopes = [];
text(0.4, 100,'(a)','fontsize',14,'fontweight','normal');
set(gca,'fontsize',12);

subplot(1,2,2);
y = BinLoc .* mean(Count(:, times(ii)-10:times(ii)+10),2);
[slope, constant] = compute_slope(BinLoc, y);
slopes = [slopes;slope];
loglog(BinLoc, y, 'b');
hold on;
loglog(BinLoc, exp(log(BinLoc)*slope + constant),'k');

load('../data/QG/hist_randph.mat');
y = BinLoc .* mean(Count(:, times(ii)-10:times(ii)+10),2);
loglog(BinLoc, BinLoc .* mean(Count(:, times(ii)-10:times(ii)+10),2),'r');
[slope, constant] = compute_slope(BinLoc, y);
slopes = [slopes;slope];
%loglog(BinLoc, exp(log(BinLoc)*slope + constant));
grid on;
legend('QG flow', ['slope: ' num2str(slope, '%.2f')], 'synthetic flow');
tl = title('$ft=1000$');
set(tl,'interpreter','latex','fontsize',14);
ylim([1e-4, 1e2]);
set(gca,'fontsize',12);
text(0.4, 100,'(b)','fontsize',14,'fontweight','normal');

print -f1 -r300 -depsc ../../../source/figures/fig2







