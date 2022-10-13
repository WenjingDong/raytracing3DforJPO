function [] = time_stepping(x, y, z, kx, ky, kz, t, dt, nt, savefreq, nsave, alpha, f, N, psihat, psihatabs, rna, rnb, ikxg, ikyg, ikzg)
[nx, ny, nz] = size(psihat);

ughat = -ikyg .* psihat;
vghat =  ikxg .* psihat;
uxghat = ikxg .* ughat;
uyghat = ikyg .* ughat;
uzghat = ikzg .* ughat;
vxghat = ikxg .* vghat;
vzghat = ikzg .* vghat;

for ii=1:nt
    
    [u1, v1, ux1, uy1, uz1, vx1, vz1] = rhs(ikxg, ikyg, ikzg, ...
        ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, x, y, z);
    vy1 = - ux1;
    cg1 = dispersion3D(kx, ky, kz, f, N);
    omega = sqrt(N^2 * (kx.^2 + ky.^2) + f^2 .* kz.^2)./sqrt(kx.^2 + ky.^2 + kz.^2) ...
                 + u1 .* kx + v1 .* ky; 
             
    xint = x + (u1 + cg1(1,:)) * dt;
    yint = y + (v1 + cg1(2,:)) * dt;
    zint = z + (0  + cg1(3,:)) * dt;
    kxint = kx - (ux1 .* kx + vx1 .* ky) * dt;
    kyint = ky - (uy1 .* kx + vy1 .* ky) * dt;  
    kzint = kz - (uz1 .* kx + vz1 .* ky) * dt;
    
    % update velocity at t+dt
   [ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, rna, rnb, ~] = velocity_ou3D(rna,rnb, alpha, psihatabs, dt, nx, ...
    ny, nz, ikxg, ikyg, ikzg);

    [u2, v2, ux2, uy2, uz2, vx2, vz2] = rhs(ikxg, ikyg, ikzg, ...
        ughat, vghat, uxghat, uyghat, uzghat, vxghat, vzghat, xint, yint, zint);
    vy2 = - ux2;
    cg2 = dispersion3D(kxint, kyint, kzint, f, N);
    
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
                       
    t = t + dt;  
    x = x2; y = y2; z = z2;
    kx = kx2; ky = ky2; kz = kz2;
    % save data
    if mod(ii,savefreq)==0
        savedata(ii/savefreq+1, nsave, t, [kx', ky', kz', x', y', z', omega'],1);       
        %if ii/savefreq + 1 == nsave
	    %   save('output.mat','t','rna','rnb','kx','ky','kz','x','y','z');
       %end	       
    end  
    
end


