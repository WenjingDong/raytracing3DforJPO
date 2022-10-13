function [kz] = reflectBC(z, kz, cgz, dt, Lz)
zint = kz + cgz/2 * dt;
index = (zint + Lz < 0) | (zint > 0);
kz(index) = -kz(index);





