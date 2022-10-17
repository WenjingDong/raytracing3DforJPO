clear;
rescale = 0;% rescale vertical length scale or not.... For the simulations in the paper, rescale=0. rescale happens in the main program.
load('psiN100.mat');
[nx, ny, nz] = size(psihat);
N = 100; f = 1;
zetahat = vxghat - uyghat;
Lx = 2 * pi;
Ly = Lx;
if rescale==1
        Lz = Lx * 2 * 100/1000; % 100 is from the rescale factor
else
	Lz = 2*Lx/1000;
end
[kxg, kyg, kzg] = ndgrid([0:nx/2-1 -nx/2:1:-1], [0:ny/2-1 -ny/2:1:-1], [0:nz/2-1 -nz/2:1:-1]);
kxg = kxg * 2*pi/Lx; kyg = kyg*2*pi/Ly; kzg = kzg*2*pi/Lz;
ikxg = 1i * kxg; ikyg = 1i * kyg; ikzg = 1i * kzg;

randph = 2 * pi * rand([nx, ny, nz]);
randph(nx:-1:nx/2+2, 1, 1:nz/2) = - randph(2:nx/2, 1, 1:nz/2);
randph(1, ny:-1:ny/2+2, :) = -randph(1, 2:ny/2, :);
randph(2:nx, ny/2+2:ny, :) = -randph(nx:-1:2, ny/2:-1:2,:);
randph(1:nx, 1:ny, nz:-1:nz/2+2) =  randph(1:nx, 1:ny, 2:nz/2); 

psihat_randph = generate_QGrandphhelper(randph, psihat);
psi_randph = ifftn(psihat_randph,'symmetric') *nx*ny*nz;
%psihat_randph = fftn(psi_randph)/(nx*ny*nz);

zetahat_randph = (ikxg.^2 + ikyg.^2) .* psihat_randph;
zeta_randph = ifftn(zetahat_randph,'symmetric')*nx*ny*nz;
rms_zeta_randph = sqrt(mean(mean(mean( zeta_randph.^2 ,1), 2), 3));
zeta = ifftn(zetahat,'symmetric')*nx*ny*nz;
rms_zeta = sqrt(mean(mean(mean( zeta.^2 ,1), 2), 3));

disp(rms_zeta_randph);
disp(rms_zeta);

% save data
psihat = psihat_randph;
ughat = -ikyg .* psihat;
vghat =  ikxg .* psihat;
uxghat = ikxg .* ughat;
uyghat = ikyg .* ughat;
vxghat = ikxg .* vghat;
vyghat = -uxghat;
uzghat = ikzg .* ughat;
vzghat = ikzg .* vghat;

u = real(ifftn(ughat))*nx*ny*nz;
v = real(ifftn(vghat))*nx*ny*nz;
ux = real(ifftn(uxghat))*nx*ny*nz;
uy = real(ifftn(uyghat))*nx*ny*nz;
uz = N/f*real(ifftn(uzghat))*nx*ny*nz;
vx = real(ifftn(vxghat))*nx*ny*nz;
vy = -ux;
vz = N/f*real(ifftn(vzghat))*nx*ny*nz;

%save('psi_randph.mat', 'psihat', 'ughat', 'vghat', 'uxghat', 'uyghat', 'uzghat', 'vxghat', 'vyghat', 'vzghat');
save('../../postprocess/data/velocity_eigen_analysis/velocity_randph.mat', 'u','v','ux','uy','uz','vx','vy','vz');







