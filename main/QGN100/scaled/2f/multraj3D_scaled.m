% Multiple trajectories with 3D velocity
set(0, 'DefaultLineLineWidth', 1);
clear;
close;
rng('shuffle');

restart = 0;
f = 1;
N = 100*f;
omega0 = 2*f;
kh0 = 50;
costh = sqrt(N^2 - omega0^2)/sqrt(N^2 - f^2); 
K0 = kh0/sqrt(1 - costh^2);
m0 = K0 * costh * f/N; % non-dimensional version
%m0_ndim = m0 ;
%[cg0] = dispersion3D(kh0,0,m0,f,N);
[cg0] = dispersion3D_ndim(kh0, 0, m0, f, N);
disp(sqrt(sum(cg0.^2)));

file = dir('psiN*.mat');
load([file.name]);
[nx, ny, nz] = size(ughat);
ns = 4;
Lx = 2*pi;  
Ly = Lx;
Lz = 2*Lx/1000 * N/f; 

% Generate wavenumber grid
[kxg, kyg, kzg] = ndgrid([0:nx/2-1 -nx/2:1:-1],[0:ny/2-1 -ny/2:1:-1],[0:nz/2-1 -nz/2:1:-1]);
kxg = kxg *2*pi/Lx; kyg = kyg *2*pi/Lx; kzg = kzg*2*pi/Lz;
ikxg = 1i * kxg; ikyg = 1i * kyg; ikzg = 1i * kzg;


u = ifftn(ughat,'symmetric')*nx*ny*nz;
v = ifftn(vghat,'symmetric')*nx*ny*nz;
disp('mean velocity magnitude is: ');
U = sqrt(mean(mean(mean(u.^2 + v.^2))));
disp(sqrt(mean(mean(mean(u.^2 + v.^2)))));

dt = 0.8*0.125*0.00125/max([max(abs(cg0)), U]);
nt = floor(1200/f/dt)+1;
savefreq = floor(1/f/dt);
nsave = floor(nt/savefreq) + 1;

[x,y,z,kx,ky,kz,t] = initialize(restart, ns, kh0, m0, Lx, Lz, nsave);

%[ughat, vghat, uxghat, uyghat, uzghat, vxghat, vyghat, vzghat, rna, rnb]...
%    = velocity_ou3D(rna, rnb, alpha, ceta, dt, nx, ny, nz, ikxg, ikyg, ikzg);
% check velocity

tic
for ii=1:nt
    
    [u1, v1, ux1, uy1, uz1, vx1, vz1] = rhs(ikxg, ikyg, ikzg, ...
        ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, x, y, z, Lz);
    vy1 = - ux1;
    cg1 = dispersion3D_ndim(kx, ky, kz, f, N);
    omega = N * sqrt( kx.^2 + ky.^2 + kz.^2)./sqrt(kx.^2 + ky.^2 + N^2/f^2 * kz.^2) ...
                 + u1 .* kx + v1 .* ky; 
             
    xint = x + (u1 + cg1(1,:)) * dt;
    yint = y + (v1 + cg1(2,:)) * dt;
    zint = z + (0  + cg1(3,:)) * dt;
    kxint = kx - (ux1 .* kx + vx1 .* ky) * dt;
    kyint = ky - (uy1 .* kx + vy1 .* ky) * dt;  
    kzint = kz - (uz1 .* kx + vz1 .* ky) * dt;
    
    [u2, v2, ux2, uy2, uz2, vx2, vz2] = rhs(ikxg, ikyg, ikzg, ...
        ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, xint, yint, zint, Lz);
    vy2 = - ux2;
    cg2 = dispersion3D_ndim(kxint, kyint, kzint, f, N);
    
    % stepping
    x2 = x + (u1 + cg1(1,:) + u2 + cg2(1,:))/2 * dt;
    y2 = y + (v1 + cg1(2,:) + v2 + cg2(2,:))/2 * dt;
    z2 = z + (cg1(3,:) + cg2(3,:))/2*dt;
    kx2 = kx - ( ux1 .* kx    + vx1 .* ky  ...
               + ux2 .* kxint + vx2 .* kyint) * dt/2;                       
    ky2 = ky - ( uy1 .* kx    + vy1 .* ky ...
               + uy2 .* kxint + vy2 .* kyint) * dt/2;                     
    kz2 = kz - ( uz1 .* kx    + vz1 .* ky ...
               + uz2 .* kxint + vz2 .* kyint) * dt/2;
                  
    % reflect boundary
    [z2, kz2] = reflectBC(z2, kz2, Lz);
    t = t + dt;  
    x = x2; y = y2; z = z2;
    kx = kx2; ky = ky2; kz = kz2;
    % save data
    if mod(ii,savefreq)==0
        savedata(ii/savefreq+1, nsave, t, [kx', ky', kz', x', y', z', omega'],1);       
        if ii/savefreq + 1 == nsave
	       save('output.mat','t','kx','ky','kz','x','y','z');
       end	       
    end  
    
end
toc

function [x,y,z,kx,ky,kz,t] = initialize(restart, ns, kh0, m0, Lx, Lz, nsave)
     if restart == 1
	    load('output.mat','t','kx','ky','kz','x','y','z');
     else
        t = 0;
        %[xloc, yloc] = ndgrid([0:sqrt(ns)-1],[0:sqrt(ns)-1]);
        %xloc = xloc/sqrt(ns)*Lx; yloc = yloc/sqrt(ns)*Lx;
	x = Lx * rand(1, ns); 
        y = Lx * rand(1, ns);
        z = -Lz/2 * rand(1,ns);

        randph =  2 * pi * rand(1,ns);
        kx = ones(1,ns) * kh0 .* cos(randph);
        ky = ones(1,ns) * kh0 .* sin(randph);
        kz = [ones(1,ns/2), -ones(1,ns/2)] * m0;

        savedata(1, nsave, t, [kx', ky', kz', x', y', z', zeros(1,ns)'],0);
      	%rna =  normrnd(0,1,[nx,ny, nz/2]);
        %rnb =  normrnd(0,1,[nx,ny, nz/2]);
      end
end 

function [u, v, ux, uy, uz, vx, vz] = rhs(ikxg, ikyg, ikzg, ...
    ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, x, y, z, Lz)

u = zeros(1,numel(x));
v = zeros(1,numel(x));
ux = zeros(1,numel(x));
uy = zeros(1,numel(x));
uz = zeros(1,numel(x));
vx = zeros(1,numel(x));
vz = zeros(1,numel(x));

for jj = 1:numel(x)
        exploc =  exp(ikxg .* x(jj) + ikyg .* y(jj) + ikzg * (z(jj) - Lz/2));
        u(jj) = real( sum( sum( sum( ughat .* exploc, 1 ),2 ), 3 ) );
        v(jj) = real( sum( sum( sum( vghat .* exploc, 1 ),2 ), 3 ) ); 
        ux(jj) = real( sum( sum( sum( uxghat .* exploc, 1 ),2 ),3 ) ); %u_x
        uy(jj) = real( sum( sum( sum( uyghat .* exploc, 1 ),2 ),3 ) ); %u_y
        uz(jj) = real( sum( sum( sum( uzghat .* exploc, 1 ),2 ),3 ) ); %u_z
        vx(jj) = real( sum( sum( sum( vxghat .* exploc, 1 ),2 ),3 ) );
        vz(jj) = real( sum( sum( sum( vzghat .* exploc, 1 ),2 ),3 ) ); %v_z
end

end

function [] = savedata(fnum, nfile, t, kxvec, flag)
     
     disp(['write data ' num2str(fnum) ' of ' num2str(nfile)]);
     if flag==0
         fileID = fopen('output.bin','w');
     else
         fileID = fopen('output.bin','a');
     end
     fwrite(fileID, kxvec, 'double');
     fclose(fileID);
     
     if flag==0
         fileID = fopen('time.bin','w');
     else
         fileID = fopen('time.bin','a');
     end
     fwrite(fileID, t, 'double');
     fclose(fileID);
end
function [z, m] = reflectBC(z, m, Lz)
     % treatment of the top boundary 
     index = (z>0);
     z(index) = -z(index);
     m(index) = -m(index);
     % treatment of the bottome boundary
     index = (z<-Lz/2);
     z(index) = -Lz - z(index);
     m(index) = -m(index);
end

function [cg] = dispersion3D_ndim(kx, ky, kz, f, N)

% compute intrinstic group velocity for Boussinesq
omegar = N * sqrt( kx.^2 + ky.^2 + kz.^2)./sqrt(kx.^2 + ky.^2 + N^2/f^2 * kz.^2);
cg(1,:) = N^2/f^2 * (N^2 - f^2)./omegar .* kx .* kz.^2 ./(kx.^2 + ky.^2 + N^2/f^2 * kz.^2).^2;
cg(2,:) = N^2/f^2 * (N^2 - f^2)./omegar .* ky .* kz.^2 ./(kx.^2 + ky.^2 + N^2/f^2 * kz.^2).^2;
cg(3,:) = N^2/f^2 * (f^2 - N^2)./omegar .* kz .* (kx.^2 + ky.^2) ./(kx.^2 + ky.^2 + N^2/f^2 * kz.^2).^2;
% determinstic velocity
%cgx = 0.01; cgy = 1;
end
