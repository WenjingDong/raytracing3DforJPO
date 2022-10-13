function [BinLoc, pdfValues, omegarFilt, DiscMatFilt] = eigen_hist(f, N, ux, uy, uz, vx, vy, vz)
% compute histogram using the eigenvector analysis
[nx, ny, nz ] = size(ux);
omegar = zeros(nx, ny, nz);
DiscMat = ux.^2 + vx .* uy; 

for ii=1:nx
	for jj=1:ny
	    for mm = 1:nz
                S = -[ux(ii,jj,mm),uy(ii,jj,mm),uz(ii,jj,mm);vx(ii,jj,mm),vy(ii,jj,mm),vz(ii,jj,mm);0,0,0]';
                [V,D] = eig(S);        
                ind = find(abs(imag(diag(D)))<1e-6 & real(diag(D))>1e-10);
                if length(ind)>1
                    disp('found more than 2 positive  eigenvalues');
                    return;  
                end
                if length(ind)>0 
	            omegar(ii,jj,mm) = sqrt(N^2 * V(1,ind)^2 + N^2 * V(2,ind)^2 + f^2 * V(3,ind)^2)/...
	          sqrt(V(1,ind)^2 + V(2,ind)^2 + V(3,ind)^2);
                 end
            end
       end
end

BinWidth = 0.02;
BinEdges = [log10(f) - BinWidth/2 : BinWidth : log10(N) + BinWidth/2];
BinLoc = 10.^( (BinEdges(1:end-1) + BinEdges(2:end))/2 );
BinEdges = 10.^BinEdges;
Width = BinEdges(2:end) - BinEdges(1:end-1);


index = abs(imag(omegar))<1e-20 & omegar > 0.5;
omegarFilt = real(omegar(index));

pdfValues = histcounts(omegarFilt,BinEdges);
pdfValues = pdfValues./Width/sum(pdfValues);

DiscMatFilt = DiscMat(DiscMat>0);

end
