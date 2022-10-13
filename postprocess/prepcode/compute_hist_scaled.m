clear;
%close all;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesTitleFontWeight','normal');
prompt1 = 'what is the value of f\n';
f = input(prompt1); 
prompt2 = 'what is the value of N\n';
N = input(prompt2); 

BinWidth = 0.1*f;%*N/20;
logbin = 1;

ns = 4;
ndata = 1000;

dirname = '/scratch/wd583/raytracing3D/arrayjobs/QGN100/scaled/';
folders = dir('/scratch/wd583/raytracing3D/arrayjobs/QGN100/scaled/randph_run*');
nrun = length(folders);
x = zeros(nrun,ns,ndata);
y = zeros(nrun,ns,ndata);
z = zeros(nrun,ns,ndata);
k = zeros(nrun,ns,ndata);
l = zeros(nrun,ns,ndata);
m = zeros(nrun,ns,ndata);
omega = zeros(nrun,ns,ndata-1);

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
    x(ii,:,:) = squeeze(output(:, 4, 1:ndata));
    y(ii,:,:) = squeeze(output(:, 5, 1:ndata));
    z(ii,:,:) = squeeze(output(:, 6, 1:ndata));

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
x_2d = reshape(x,[nrun*ns,ndata]);
y_2d = reshape(y,[nrun*ns,ndata]);
z_2d = reshape(z,[nrun*ns,ndata]);
k_2d = reshape(k,[nrun*ns,ndata]);
l_2d = reshape(l,[nrun*ns,ndata]);
m_2d = reshape(m,[nrun*ns,ndata]);
K_2d = sqrt(k_2d.^2 + l_2d.^2 + m_2d.^2);
kh_2d= sqrt(k_2d.^2 + l_2d.^2);
omegar_2d = N * sqrt( k_2d.^2 + l_2d.^2 + m_2d.^2)./sqrt(k_2d.^2 + l_2d.^2 + N^2/f^2 * m_2d.^2);

kh_2d = sqrt(k_2d.^2 + l_2d.^2);
cg = (N^2 - f^2).* kh_2d .* abs(m_2d)./(K_2d.^3 .* omegar_2d);

% compute time series of histcounts of intrinsic frequency
%Edges = [f-BinWidth/2 : BinWidth : N + BinWidth/2]';
%BinLoc = (Edges(1:end-1) + Edges(2:end))/2;
Edges = [log10(f)-0.01 : 0.02 : log10(N)+0.01]';
BinLoc = 10.^((Edges(1:end-1) + Edges(2:end))/2);
Edges = 10.^Edges;
Count = zeros(length(Edges)-1,ndata);
Width = Edges(2:end) - Edges(1:end-1);

for kk=1:ndata
	Count(:,kk) = histcounts(omegar_2d(:,kk),Edges);
        Count(:,kk) = Count(:,kk)./Width/(nrun*ns);
end
[BinLoc_2d,t_2d] = ndgrid(BinLoc,t);

times = [100, 250, 600, 800, 1000];
BinLoc = BinLoc/f;

figure('unit','centimeter','position',[15 35 27 7]);
subplot(1,2,1);
loglog(BinLoc, BinLoc.*Count(:,1),'*');
hold on;
loglog(BinLoc,  BinLoc.*mean(Count(:,times(1)-10:times(1)+10),2)); 
hold on;
loglog(BinLoc,  BinLoc.*mean(Count(:,times(2)-10:times(2)+10),2)); 
loglog(BinLoc,  BinLoc.*mean(Count(:,times(3)-10:times(3)+10),2)); 
loglog(BinLoc,  BinLoc.*mean(Count(:,times(4)-10:times(4)+10),2)); 
loglog(BinLoc,  BinLoc.*mean(Count(:,times(5)-20:times(5)),2)); 
grid on;
xlabel('$\omega_r/f$','interpreter','latex','fontsize',12);
ylabel('$e(\omega_r)$','interpreter','latex','fontsize',12);
tl = title('Time evolution of energy density $e(\omega_r)$');
set(tl,'interpreter','latex');
labels = {'$ft$ = 0 ',...
['$ft$ = ' num2str(f*t(times(1)),'%.0f')],...
['$ft$ = ' num2str(f*t(times(2)),'%.0f')],...
['$ft$ = ' num2str(round(f*t(times(3))))],...
['$ft$ = ' num2str(round(f*t(times(4))))],...
['$ft$ = ' num2str(round(f*t(times(5))))]};
lg = legend(labels);
set(lg,'interpreter','latex');
xlim([1,N/f]);
set(gca,'xtick',[1,2,10 20,40,N/f]);
%ylim([1e-3,2/BinWidth*f]);
set(gca,'fontsize',10);

subplot(1,2,2);
loglog(BinLoc, BinLoc.*mean(Count(:,times(5)-20:times(5)),2)); 
hold on;
grid on;
%loglog(BinLoc, BinLoc.*mean(Count(:,ndata/5-25:ndata/5+25),2)); 
y = log(BinLoc.*mean(Count(:,ndata-20:ndata),2)); 
index1 = (BinLoc > 2) & (BinLoc < min([N/f*0.5,40]));
index2 = ~(isnan(y) | isinf(y)); 
index = index1 & index2;
x = log(BinLoc);
p = polyfit(x(index), y(index), 1);
ypol = polyval(p,x);
loglog(BinLoc(index),exp(ypol(index)),'k');

y1 = log(BinLoc.*mean(Count(:,ndata/5-10:ndata/5+10),2));
p1 = polyfit(x(index),y1(index),1);
ypol1 = polyval(p1,x);
%loglog(BinLoc(index),exp(ypol1(index)),'r');

title(['slope: ' num2str(p(1),'%.2f')]);
%ylim([1e-3,1e1]);
xlim([1,N/f]);
set(gca,'xtick',[1,2,10 20,40,N/f]);
ylabel('Time averaged $e(\omega_r)$','interpreter','latex','fontsize',12);


%save('data/QG/hist_randph.mat', 't', 'BinLoc', 'Count');

