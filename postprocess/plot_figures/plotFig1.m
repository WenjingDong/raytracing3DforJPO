clear;
close all;
Lx = 2 * pi; Ly = Lx;

load('../data/QGflows/psiN100.mat');
[nx, ny, nz] = size(psihat);
zetahat = vxghat - uyghat;
zeta = ifftn(zetahat, 'symmetric') * nx * ny * nz;
zeta_qg_mid = ( zeta(:,:,(nz+2)/4) + zeta(:,:,(nz+2)/4+1) )/2;
[x, y] = ndgrid([0:nx-1]/nx*Lx, [0:ny-1]/ny*Ly);
x_ld = x * 5; y_ld = y * 5; 
disp(rms(zeta(:)));
clear zetahat;

load('../data/QGflows/psi_randph.mat');
zetahat = vxghat - uyghat;
zeta = ifftn(zetahat, 'symmetric') * nx * ny * nz;
zeta_rand_mid = ( zeta(:,:,(nz+2)/4) + zeta(:,:,(nz+2)/4+1) )/2;
disp(rms(zeta(:)));

figure('unit', 'centimeter', 'position', [10, 10, 27, 10]);
subplot(1,2,1);
pcolor(x_ld, y_ld, zeta_qg_mid');
shading interp;colorbar;
tl = title('QG flow'); 
set(tl, 'interpreter', 'latex', 'fontsize', 14);
MAX = max([max(max(abs(zeta_qg_mid))), max(max(abs(zeta_rand_mid)))]);
caxis([-MAX, MAX]);
axis image;
xlabel('$x/L_d$', 'interpreter', 'latex', 'fontsize', 14);
ylabel('$y/L_d$', 'interpreter', 'latex', 'fontsize', 14);
text(-5, 32, '(a)', 'fontsize', 14);
colormap(redblue(100));
set(gca,'fontsize',12);

subplot(1,2,2);
pcolor(x_ld, y_ld, zeta_rand_mid');
shading interp;colorbar;
tl = title('Synthetic flow'); 
set(tl, 'interpreter', 'latex', 'fontsize', 14);
MAX = max(max(abs(zeta_rand_mid)));
caxis([-MAX, MAX]);
axis image;
text(-5, 32, '(b)','fontsize',14);
colormap(colormap);
set(gca,'fontsize',12);

print -f1 -depsc -r300 ../../source/figures/qg_randph_zeta
