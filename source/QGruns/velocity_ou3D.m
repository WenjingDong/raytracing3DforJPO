function [ughat, vghat, uxghat, uyghat, uzghat, vxghat, vyghat, vzghat, rna1, rnb1] = velocity_ou3D(rna,rnb, alpha, ceta, dt, NX, ...
    NY, NZ, ikxg, ikyg, ikzg)

% rn: a_k, b_k at previous step, rn1: next step
rna1 = sqrt(1 - exp(-2*alpha*dt)) .* normrnd(0,1,[NX,NY,NZ/2]);
rna1 = rna1 +  rna * exp(-alpha * dt);

rnb1 = sqrt(1 - exp(-2*alpha*dt)) .* normrnd(0,1,[NX,NY, NZ/2]);
rnb1 = rnb1 +  rnb * exp(-alpha * dt);


psihat =  ceta .* (rna1 + 1i * rnb1);

psihat(:,NY/2+1,1:NZ/2) = 0;
psihat(NX/2+1,:,1:NZ/2) = 0;
psihat(:,:,NZ/2+1) = 0;

psihat(:,:,1) = sympad_(psihat(:,1:NY/2,1),0);

psihat(NX:-1:2, NY:-1:2, NZ:-1:NZ/2+2) = conj(psihat(2:end,2:end,2:NZ/2));

psihat(1, NY:-1:2, NZ:-1:NZ/2+2) = conj(psihat(1, 2:NY, 2:NZ/2));

psihat(NX:-1:2, 1, NZ:-1:NZ/2+2) = conj(psihat(2:NX, 1, 2:NZ/2));

psihat(1, 1, NZ:-1:NZ/2+1) = conj(psihat(1, 1, 1:NZ/2));

%psi = ifft2(psihat,'symmetric')*(NX^2);

ughat = -ikyg .* psihat;
vghat =  ikxg .* psihat;
uxghat = ikxg .* ughat;
uyghat = ikyg .* ughat;
uzghat = ikzg .* ughat;
vxghat = ikxg .* vghat;
vyghat = - uxghat;
vzghat = ikzg .* vghat;
 
end