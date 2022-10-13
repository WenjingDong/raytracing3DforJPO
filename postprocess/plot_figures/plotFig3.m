clear;
close all;
set(0,'DefaultLineLineWidth',2);
files = dir('data/velocity_eigen_analysis/velocity_*.mat');
fig3 = figure('unit','centimeter','position', [10, 10, 12, 8] );
for ii=1:2
	load(['data/velocity_eigen_analysis/' files(ii).name]);
	D = ux.^2 + vx.* uy;
        figure(ii);	
	h = histogram(D(:));

        figure(3);
	BinLoc = (h.BinEdges(1:end-1) + h.BinEdges(2:end))/2;
        plot(BinLoc, h.Values);
        hold on;
end
xlim([-0.02, 0.02]);
grid on;
legend('QG flow', 'Synthetic flow');
set(gca, 'fontsize',12);
xlabel('$D/f$','interpreter','latex','fontsize',14);
ylabel('p.d.f. of $D$', 'interpreter','latex', 'fontsize', 14);
close([1, 2]);
print -f3 -depsc -r300 ../source/figures/hist_D
