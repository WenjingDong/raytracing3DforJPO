function [x,y,z,kx,ky,kz,t] = initialize(restart, ns, kh0, m0, Lx, Lz, nx, ny, nz, nsave)
     if restart == 1
	    load('output.mat','rna','rnb','t','kx','ky','kz','x','y','z');
     else
        t = 0;
        rng('shuffle'); 
        x = Lx * rand(1, ns); y = Lx * rand(1, ns);
	    z = -Lz * rand(1, ns);

        randph = 2 * pi * rand(1,ns);%rand(1,ns);
        kx = ones(1,ns) * kh0 .* cos(randph);
        ky = ones(1,ns) * kh0 .* sin(randph);
        kz = [ones(1,ns/2), -ones(1,ns/2)] * m0;

        savedata(1, nsave, t, [kx', ky', kz', x', y', z', zeros(1,ns)'],0);
      end
end 
