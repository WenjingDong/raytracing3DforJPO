function [u, v, ux, uy, uz, vx, vz] = rhs(ikxg, ikyg, ikzg, ...
    ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, x, y, z)
% evalute velocity of rays for given locations (x, y, z),
% and the fourier series of velocity field 
u = zeros(1,numel(x));
v = zeros(1,numel(x));
ux = zeros(1,numel(x));
uy = zeros(1,numel(x));
uz = zeros(1,numel(x));
vx = zeros(1,numel(x));
vz = zeros(1,numel(x));

for jj = 1:numel(x)
        exploc =  exp(ikxg .* x(jj) + ikyg .* y(jj) + ikzg * z(jj));
        u(jj) = real( sum( sum( sum( ughat .* exploc, 1 ),2 ), 3 ) );
        v(jj) = real( sum( sum( sum( vghat .* exploc, 1 ),2 ), 3 ) ); 
        ux(jj) = real( sum( sum( sum( uxghat .* exploc, 1 ),2 ),3 ) ); %u_x
        uy(jj) = real( sum( sum( sum( uyghat .* exploc, 1 ),2 ),3 ) ); %u_y
        uz(jj) = real( sum( sum( sum( uzghat .* exploc, 1 ),2 ),3 ) ); %u_z
        vx(jj) = real( sum( sum( sum( vxghat .* exploc, 1 ),2 ),3 ) );
        vz(jj) = real( sum( sum( sum( vzghat .* exploc, 1 ),2 ),3 ) ); %v_z
end

end
