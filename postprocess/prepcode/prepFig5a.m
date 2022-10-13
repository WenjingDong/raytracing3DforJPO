% This script computes average intrinsic frequency used for figure 5 (a)
clear;
close all;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesTitleFontWeight','normal');
prompt1 = 'what is the value of f\n';
f = input(prompt1); 
prompt2 = 'what is the value of N\n';
N = input(prompt2); 

BinWidth = 0.1*f;%*N/20;
ns = 100;
ndata = 5600;

dirname = '/scratch/wd583/raytracing3D/arrayjobs/charney/N100/5f/';
folders = dir('/scratch/wd583/raytracing3D/arrayjobs/charney/N100/5f/5f_run*');
nrun = length(folders);
k = zeros(nrun,ns,ndata);
l = zeros(nrun,ns,ndata);
m = zeros(nrun,ns,ndata);
omega = zeros(nrun,ns,ndata-1);

for ii=1:nrun
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
dt = t(2) - t(1);
K  = sqrt(k.^2 + l.^2 + m.^2);
kh = sqrt(k.^2 + l.^2);
omegar = sqrt(N^2 * (k.^2 + l.^2) + f^2 * m.^2)./K;

omegar_2d = reshape(omegar,[nrun*ns,ndata]);
omegar_std = std(omegar_2d,0,1);
omegar_avg = squeeze(mean(omegar_2d, 1));
plot(t, omegar_std);
save('data/charney/omegaavg_init_Ro01f5.mat', 't', 'omegar_avg');
