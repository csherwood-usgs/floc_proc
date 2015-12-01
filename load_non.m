% load_non.m - Load the inactive classes in runs > 66

snn = zeros( NNN, nz, nt);
for n=1:NNN
    ncname = sprintf('sand_%02d',n)
    snn(n,:,:)=squeeze(ncread(url,ncname,[i j 1 1],[1 1 Inf Inf]));
end
snns = squeeze(sum(snn));
sumsnns = sum( snns.*dzw);

% total conc
figure(5); clf
pcolorjw( s2d*tz, h+z_w, log10(snns+eps))
colorbar
title('Non-depositing Sand log_{10} Total Concentration')

% fraction-weighted ws
wst = squeeze(sum(repmat(ws(NCS+1:NCS+NNN),1,nz,nt).*snn)./sum(snn));
figure(6); clf
pcolorjw( s2d*tz, h+z_w, 1e3*wst)
colorbar
title('Non-depositing Sand Settling Velocity (mm/s)')

% fraction-weighted size
figure(7); clf
fdiamt = squeeze(sum(repmat(fdiam(NCS+1:NCS+NNN),1,nz,nt).*snn)./sum(snn));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamt)
colorbar
title('Non-depositing Sand Diameter (\mum)')

% acoustic response
Dfv = fdiam(NCS+1:NCS+NNN);
rhofv = rhos(NCS+1:NCS+NNN);
vnn = zeros(nz,nt);
for jj=1:nt
    for ii=1:nz
        mv = squeeze(snn(:,ii,jj));
        vnn(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,3e6);
    end
end
figure(8)
pcolorjw( s2d*tz, h+z_w, vnn)
colorbar
title('Non-depositing Acoustic Response')
xlabel('Days')
ylabel('Elevation (m)')
xlim([0 32.8])
caxis([ 0 .3])
ylim([0 3])