function [psihat] = generate_QGrandphhelper(randph, ceta)

psihat =  ceta .* exp(1i * randph);
[NX, NY, NZ] = size(psihat);

psihat(:,NY/2+1,1:NZ/2) = 0;
psihat(NX/2+1,:,1:NZ/2) = 0;
psihat(:,:,NZ/2+1) = 0;
psihat(:,:,1) = sympad_(psihat(:,1:NY/2,1),0);
psihat(NX:-1:2, NY:-1:2, NZ:-1:NZ/2+2) = conj(psihat(2:end,2:end,2:NZ/2));
psihat(1, NY:-1:2, NZ:-1:NZ/2+2) = conj(psihat(1, 2:NY, 2:NZ/2));
psihat(NX:-1:2, 1, NZ:-1:NZ/2+2) = conj(psihat(2:NX, 1, 2:NZ/2));
psihat(1, 1, NZ:-1:NZ/2+1) = conj(psihat(1, 1, 1:NZ/2));

end
