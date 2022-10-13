clear;
%close all;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesTitleFontWeight','normal');
prompt1 = 'what is the value of f\n';
f = input(prompt1); 
prompt2 = 'what is the value of N\n';
N = input(prompt2); 
prompt3 = 'logbin?\n';
logbin = input(prompt3);

BinWidth = 0.1*f;%*N/20;
ns = 100;
ndata = 1200;

dirname = '/scratch/wd583/raytracing3D/arrayjobs/charney_gaussian/N20_unsteady/Ro01/';
folders = dir('/scratch/wd583/raytracing3D/arrayjobs/charney_gaussian/N20_unsteady/Ro01/2f_run*');
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
omegar_2d = reshape(omegar,[nrun*ns, ndata]);
omega_2d = reshape(omega, [nrun*ns, ndata-1]);

% compute time series of histcounts of intrinsic frequency
if logbin==0
   Edges = [f-BinWidth/2 : BinWidth : N + BinWidth/2]';
   BinLoc = (Edges(1:end-1) + Edges(2:end))/2;
else
   Edges = [log10(f)-0.01 : 0.02 : log10(N)+0.01]';
   BinLoc = 10.^((Edges(1:end-1) + Edges(2:end))/2);
   Edges = 10.^Edges;
end

Count = zeros(length(Edges)-1,ndata-1);
Width = Edges(2:end) - Edges(1:end-1);

for kk=1:length(Width)
        indicator = omega_2d>=Edges(kk) & omega_2d<Edges(kk+1);	
	Count(kk,:) = sum(indicator .* omegar_2d(:,2:end), 1); 
        Count(kk,:) = Count(kk,:)/Width(kk)/nrun/ns;
end

times = [ndata/10+1, ndata/5+1, ndata/5*3+1, ndata/5*4+1, ndata-10];
BinLoc = BinLoc/f;

%figure('unit','centimeter','position',[15 35 27 8]);
subplot(1,2,2);
loglog(BinLoc, BinLoc.*Count(:,1),'*');
hold on;
loglog(BinLoc,  BinLoc.*mean(Count(:,times(1)-10:times(1)),2)); 
hold on;
loglog(BinLoc,  BinLoc.*mean(Count(:,times(2)-10:times(2)),2)); 
loglog(BinLoc,  BinLoc.*mean(Count(:,times(3)-10:times(3)),2)); 
loglog(BinLoc,  BinLoc.*mean(Count(:,times(4)-10:times(4)),2)); 
loglog(BinLoc,  BinLoc.*mean(Count(:,times(5)-10:times(5)),2)); 
grid on;
xlabel('$\omega_r/f$','interpreter','latex','fontsize',12);
ylabel('$e(\omega_r)$','interpreter','latex','fontsize',12);
tl = title('Time evolution of energy density $e(\omega_r)$');
set(tl,'interpreter','latex');
labels = {'$ft$ = 0 ',...
['$ft$ = ' num2str(f*t(times(1)),'%.0f')],...
['$ft$ = ' num2str(f*t(times(2)),'%.0f')],...
['$ft$ = ' num2str(f*t(times(3)),'%.0f')],...
['$ft$ = ' num2str(f*t(times(4)),'%.0f')],...
['$ft$ = ' num2str(f*t(times(5)),'%.0f')]};
lg = legend(labels);
set(lg,'interpreter','latex');
xlim([1,N/f]);
set(gca,'xtick',[1:1:10 20:20:N/f]);
%ylim([1e-1,2/BinWidth*f]);
set(gca,'fontsize',10);

return;
subplot(1,2,2);
loglog(BinLoc, BinLoc.*mean(Count(:,ndata-20:ndata-1),2)); 
hold on;
grid on;
%loglog(BinLoc, BinLoc.*mean(Count(:,ndata/5-25:ndata/5+25),2)); 
y = log(BinLoc.*mean(Count(:,ndata-20:ndata-1),2)); 
index1 = (BinLoc > 2) & (BinLoc < min([N/f*0.5,40]));
index2 = ~(isnan(y) | isinf(y)); 
index = index1 & index2;
x = log(BinLoc);
p = polyfit(x(index), y(index), 1);
ypol = polyval(p,x);
semilogx(BinLoc(index),exp(ypol(index)),'k');

y1 = log(BinLoc.*mean(Count(:,ndata/5-10:ndata/5+10),2));
p1 = polyfit(x(index),y1(index),1);
ypol1 = polyval(p1,x);
%loglog(BinLoc(index),exp(ypol1(index)),'r');

title(['slope: ' num2str(p(1),'%.2f')]);
%ylim([1e-3,1e1]);
xlim([1,N/f]);
set(gca,'xtick',[1:1:10 20:20:N/f]);
ylabel('Time averaged $e(\omega_r)$','interpreter','latex','fontsize',12);

%save('data/charney/hist_Ro0025_f2.mat', 't', 'BinLoc', 'Count');

