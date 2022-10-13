function [ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, rna1, rnb1, psihat] = velocity_ou3D(rna,rnb, alpha, psihatabs, dt, nx, ...
    ny, nz, ikxg, ikyg, ikzg)
% This function update a streamfunction whose fourier transform are specified through OU process

% rn: a_k, b_k at previous step, rn1: next step
rna1 = sqrt(1 - exp(-2*alpha*dt)) .* normrnd(0,1,[nx,ny,nz]);
rna1 = rna1 +  rna * exp(-alpha * dt);

rnb1 = sqrt(1 - exp(-2*alpha*dt)) .* normrnd(0,1,[nx,ny,nz]);
rnb1 = rnb1 +  rnb * exp(-alpha * dt);

psihat =  psihatabs .* (rna1 + 1i * rnb1);

psihat(:,ny/2+1,:) = 0;
psihat(nx/2+1,:,:) = 0;
psihat(:,:,nz/2+1) = 0;

psihat(:,:,1) = sympad_(psihat(:,1:ny/2,1),0);

psihat(nx:-1:2, ny:-1:2, nz:-1:nz/2+2) = conj(psihat(2:end,2:end,2:nz/2));

psihat(1, ny:-1:2, nz:-1:nz/2+2) = conj(psihat(1, 2:ny, 2:nz/2));

psihat(nx:-1:2, 1, nz:-1:nz/2+2) = conj(psihat(2:nx, 1, 2:nz/2));

psihat(1, 1, nz:-1:nz/2+1) = conj(psihat(1, 1, 1:nz/2));


ughat = -ikyg .* psihat;
vghat =  ikxg .* psihat;
uxghat = ikxg .* ughat;
uyghat = ikyg .* ughat;
uzghat = ikzg .* ughat;
vxghat = ikxg .* vghat;
vzghat = ikzg .* vghat;
 
end
