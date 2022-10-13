% Multiple trajectories with 3D velocity
set(0, 'DefaultLineLineWidth', 1);
clear;
close;
rng(1);

restart = 0;
f = 1;
N = 100;
omega0 = 2*f;
kh0 = 10;
m0 = sqrt(N^2 - omega0^2)/sqrt(omega0^2 - f^2) * kh0;
[cg0] = dispersion3D(kh0,0,m0,f,N);
ns = 100;

Lx = 2*pi*2;  
Ly = Lx;
Lz = f/N * Lx; 
nx = 64;
ny = nx;
nz = nx;

% Generate wavenumber grid
[kxg, kyg, kzg] = ndgrid([0:nx/2-1 -nx/2:1:-1],[0:ny/2-1 -ny/2:1:-1],[0:nz/2-1 -nz/2:1:-1]);
kxg = kxg *2*pi/Lx; kyg = kyg *2*pi/Lx; kzg = kzg*2*pi/Lz;
ikxg = 1i * kxg; ikyg = 1i * kyg; ikzg = 1i * kzg;

% OU process parameters
n = 7; % exponent of spectra
Ro = 0.1;
alpha = 0.2;

kmin = 1; kmax = 50;% spectra
sgnf = zeros(nx, ny, nz);
Khg = sqrt(kxg.^2 + kyg.^2);
Ksg = sqrt(kxg.^2 + kyg.^2 + f^2/N^2 * kzg.^2); % kz scaled
Ksginv = 1./Ksg; Ksginv(1,1,1) = 0;
Khginv = 1./Khg; Khginv(1,1,:) = 0;
index = (Ksg>kmin)&(Ksg<kmax);
sgnf(index) = 1;

C =  sqrt(Ro^2/( sum( sum( sum( Khg.^4 .* sgnf .* Ksginv.^n ) ) ) ));
psihatabs = C * sqrt( 1/2 * Ksginv.^n ) .* sgnf;
%random numbers for a_k and b_k
rna = normrnd(0, 1, [nx, ny, nz]);
rnb = normrnd(0, 1, [nx, ny, nz]);

dt = 0.001/max([max(abs(cg0)), alpha]);
nt = floor(2000/dt)+1;
savefreq = floor(1/dt);
nsave = floor(nt/savefreq) + 1;

[ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, rna, rnb, psihat] ...
    = velocity_ou3D(rna, rnb, alpha, psihatabs, 0, nx, ny, nz, ikxg, ikyg, ikzg); 
% save velocity field
u = ifftn(ughat,'symmetric')*nx*ny*nz;
v = ifftn(vghat,'symmetric')*nx*ny*nz;  
ux = ifftn(uxghat,'symmetric')*nx*ny*nz;
vx = ifftn(vxghat,'symmetric')*nx*ny*nz;  
uy = ifftn(uyghat,'symmetric')*nx*ny*nz;
uz = ifftn(uzghat,'symmetric')*nx*ny*nz;
vz = ifftn(vzghat,'symmetric')*nx*ny*nz;

save('initial_streamfunction','psihat');
% check velocity
u = ifftn(ughat,'symmetric')*nx*ny*nz;
v = ifftn(vghat,'symmetric')*nx*ny*nz;  
U = sqrt(mean(mean(mean(u.^2+v.^2)))); 
zeta = ifftn(vxghat - uyghat, 'symmetric')*nx*ny*nz;

disp('mean velocity magnitude is: ');
disp(U);
disp('mean vorticity is: ');  
disp(sqrt(mean(mean(mean(zeta.^2)))));
disp('parameter epsilon is: ');
disp(U/sqrt(cg0(1)^2 + cg0(2)^2 + cg0(3)^2));


[x,y,z,kx,ky,kz,t] = initialize(restart, ns, kh0, m0, Lx, Lz, nx, ny, nz, nsave);

tic
time_stepping(x, y, z, kx, ky, kz, t, dt, nt, savefreq, nsave, alpha, f, N, psihat, psihatabs, rna, rnb, ikxg, ikyg, ikzg);
toc

