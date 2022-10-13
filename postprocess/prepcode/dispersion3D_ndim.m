function [cg] = dispersion3D_ndim(kx, ky, kz, f, N)
omegar = N * sqrt( kx.^2 + ky.^2 + kz.^2)./sqrt(kx.^2 + ky.^2 + N^2/f^2 * kz.^2);
cg(1,:) = N^2/f^2 * (N^2 - f^2)./omegar .* kx .* kz.^2 ./(kx.^2 + ky.^2 + N^2/f^2 * kz.^2).^2;
cg(2,:) = N^2/f^2 * (N^2 - f^2)./omegar .* ky .* kz.^2 ./(kx.^2 + ky.^2 + N^2/f^2 * kz.^2).^2;
cg(3,:) = N^2/f^2 * (f^2 - N^2)./omegar .* kz .* (kx.^2 + ky.^2) ./(kx.^2 + ky.^2 + N^2/f^2 * kz.^2).^2;
% determinstic velocity
%cgx = 0.01; cgy = 1;
end

