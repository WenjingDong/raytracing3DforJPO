function [cg] = dispersion3D(kx,ky, kz, f, N)
% compute intrinstic group velocity for inertia-gravity waves in three dimensional Boussinesq system
omegar = sqrt(N^2 * (kx.^2 + ky.^2) + f^2 * kz.^2)./sqrt(kx.^2 + ky.^2 + kz.^2);
cg(1,:) = (N^2 - f^2)./omegar .* kx .* kz.^2 ./(kx.^2 + ky.^2 + kz.^2).^2;
cg(2,:) = (N^2 - f^2)./omegar .* ky .* kz.^2 ./(kx.^2 + ky.^2 + kz.^2).^2;
cg(3,:) = (f^2 - N^2)./omegar .* kz .* (kx.^2 + ky.^2) ./(kx.^2 + ky.^2 + kz.^2).^2;
end
