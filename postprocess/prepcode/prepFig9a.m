clear;
close all;
f = 1; N = 100;
ns = 100;
ndata = 8000;
dirname = '/scratch/wd583/raytracing3D/arrayjobs/charney/N100/Ro005/';
folders = dir('/scratch/wd583/raytracing3D/arrayjobs/charney/N100/Ro005/2f_run*');
nrun = length(folders);
k = zeros(nrun,ns,ndata);
l = zeros(nrun,ns,ndata);
m = zeros(nrun,ns,ndata);

for ii=1:nrun
%    fileID = fopen([dirname '10f_run' num2str(ii) '/output.bin']);
    fileID = fopen([dirname folders(ii).name '/output.bin']);
    output = fread(fileID,'double');
    fclose(fileID);
    nt = length(output)/7/ns;

    output = reshape(output,[ns,7,nt]);
    omega(ii,:,:) = squeeze(output(:,7,2:ndata));
    k(ii,:,:) = squeeze(output(:,1,1:ndata));
    l(ii,:,:) = squeeze(output(:,2,1:ndata));
    m(ii,:,:) = squeeze(output(:,3,1:ndata));

    fileID = fopen([dirname folders(ii).name '/time.bin']);
    t = fread(fileID,'double');
    fclose(fileID);
    t(ndata+1:end) = [];
end

k_2d = reshape(k, [nrun*ns, ndata]);
l_2d = reshape(l, [nrun*ns, ndata]);
m_2d = reshape(m, [nrun*ns, ndata]);

cg_rms = zeros(ndata,1);
for ii=1:ndata
	cg = dispersion3D(k_2d(:,ii), l_2d(:,ii), m_2d(:,ii),f,N);
	cg_rms(ii) = squeeze(sqrt(mean(cg(1,:).^2+cg(2,:).^2 + cg(3,:).^2,2)));
end
plot(t, cg_rms);
save('data/charney/cg_Ro005.mat','t','cg_rms');

