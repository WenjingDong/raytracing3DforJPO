% remove higher vertical modes to satifisy no flux boundary condition 
% first, make an even extension of streamfunction in vertical direction
clear;
load('psi.mat');
[nx,ny,nz] = size(psi);
NZ = 2*nz - 2;
psi = psi/(2e6*0.05)*2*pi;
psi_ext = zeros(nx,ny,2*nz-2);
psi_ext(1:nx, 1:ny, 1) = psi(1:nx, 1:ny, nz);
psi_ext(1:nx, 1:ny, 2:nz-1) = psi(1:nx,1:ny,nz-1:-1:2);
psi_ext(1:nx, 1:ny, nz:2*nz-2) = psi(1:nx,1:ny,1:nz-1);
psi_org = psi;
clear psi;

Lx = 2*pi;
Ly = Lx; Lz = 2*Lx/1000;
[kxg, kyg, kzg] = ndgrid([0:nx/2-1 -nx/2:1:-1], [0:ny/2-1 -ny/2:1:-1], [0:NZ/2-1 -NZ/2:1:-1]);
kxg = kxg * 2*pi/Lx; kyg = kyg*2*pi/Ly; kzg = kzg*2*pi/Lz;
ikxg = 1i * kxg; ikyg = 1i * kyg; ikzg = 1i * kzg;

% remove high vertical modes
cutk = 20;
psihat = fft(psi_ext, [], 3);
psihat(:,:,cutk+1:NZ-cutk+1) = 0;
psi_filter = ifft(psihat, [], 3, 'symmetric');
psi = psi_filter;

save('psi_filter.mat','psi'); 

% scale the field to Ro=0.1
psihat = fftn(psi)/(nx*ny*NZ);
zetahat = (ikxg.^2 + ikyg.^2) .* psihat;
zeta = real(ifftn(zetahat))*nx*ny*NZ;
scale_factor = 0.1/rms(zeta(:));
psihat = scale_factor .* psihat;

uhat = -ikyg .* psihat;
vhat =  ikxg .* psihat;
uxhat = ikxg .* uhat;
uyhat = ikyg .* uhat;
uzhat = ikzg .* uhat;
vxhat = ikxg .* vhat;
vzhat = ikzg .* vhat;

save('psiN100.mat','psihat','ughat','uxghat','uyghat','uzghat','vghat','vxghat','vyghat','vzghat');

% generate a velocity field eigen vector analysis
u = real(ifftn(uhat))*nx*ny*NZ;
v = real(ifftn(vhat))*nx*ny*NZ;
ux = real(ifftn(uxhat))*nx*ny*NZ;
uy = real(ifftn(uyhat))*nx*ny*NZ;
uz = real(ifftn(uzhat))*nx*ny*NZ;
vx = real(ifftn(vxhat))*nx*ny*NZ;
vy = -ux;
vz = real(ifftn(vzhat))*nx*ny*NZ;

disp(sqrt(mean(mean(mean( u.^2+v.^2 )))));
%save('data/velocity_eigen_analysis/velocity_qg.mat','u','v','ux','uy','uz','vx','vy','vz');








