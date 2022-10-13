function [maxfreqvar, maxfreqvar_all] = maxfreqvariation(ns, ndata, N)
dirname = '/scratch/wd583/raytracing3D/arrayjobs/QGN100/scaled/';
folders = '/scratch/wd583/raytracing3D/arrayjobs/QGN100/scaled/2fddt_run*';
files = dir(folders);
nrun = length(files);
omega = zeros(nrun, ns, ndata-1);
k = zeros(nrun, ns, ndata);
l = zeros(nrun, ns, ndata);
m = zeros(nrun, ns, ndata);

for ii=1:nrun
        disp(ii);	
	fileID = fopen([dirname files(ii).name '/output.bin']);
	output = fread(fileID, 'double');
        fclose(fileID);
	nt = length(output)/7/ns;
        output = reshape(output,[ns,7,nt]);
        k(ii,:,:) = squeeze(output(:,1,1:ndata));
	l(ii,:,:) = squeeze(output(:,2,1:ndata));
	m(ii,:,:) = N*squeeze(output(:,3,1:ndata));

	omega(ii,:,:) = squeeze(output(:,7,2:ndata));
        
	fileID = fopen([dirname files(ii).name '/time.bin']);
        t = fread(fileID,'double');
        fclose(fileID);
        t(ndata+1:end) = [];

end

maxfreqvar_all = squeeze(squeeze(max(max( abs(omega - omega(:,:,1)) ,3),2)));
f = 1;
omegar = sqrt(N^2*(k.^2+l.^2) + f^2*m.^2)./sqrt(k.^2 + l.^2 + m.^2);
omegar = reshape(omegar, [nrun*ns, ndata]);

omega = reshape(omega, [nrun*ns, ndata-1]);
maxfreqvar = max(max(omega,[], 2) - min(omega, [], 2));

figure('unit','centimeter','position',[12,12,24,8]);;
subplot(1,2,1);
for ii=1:nrun*ns
	plot(t(1:end-1), (omega(ii,:) - omega(ii,1))/2);
	hold on;
end
xlabel('$ft$', 'interpreter', 'latex', 'fontsize', 14);
ylabel('$\Delta \omega_a/\omega_0$', 'interpreter', 'latex', 'fontsize', 14);
%ylim([-10,100]);
set(gca,'fontsize',12);
text(-200, 0.1, '(a)', 'fontsize', 14);

subplot(1,2,2);
for ii=1:nrun*ns
	plot(t, (omegar(ii,:) - omegar(ii,1))/2);
        disp(omegar(ii,1));
	hold on;
end
xlabel('$ft$', 'interpreter', 'latex', 'fontsize', 14);
ylabel('$\Delta \omega/\omega_0$', 'interpreter', 'latex', 'fontsize', 14);
%ylim([-10,100]);
set(gca,'fontsize',12);
text(-200, 50, '(b)', 'fontsize', 14);
%tl = title(['$\max{\Delta \omega_r}/\omega_0 = $' num2str(maxfreqvar/2,'%.2f%')]);
%set(tl, 'interpreter', 'latex', 'fontsize', 12);

%maxfreqvar_all = max(omega,[],2) - min(omega,[],2);
disp(size(maxfreqvar_all));
print -f1 -depsc -r300 ../source/figures/cons_absfreq






















