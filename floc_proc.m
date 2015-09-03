% floc_proc.m - Script to read and plot ROMS .his files
clear

cas = 66
url = sprintf('ocean_his%2d.nc',cas)
% url = 'http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/aretxabaleta/MVCO/ocean_his_44.nc'
%url = sprintf('http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/aretxabaleta/MVCO/ocean_his_%02d.nc', cas)
% Read in NST and Nbed instead of loading in huge files to infer their size
ncid = netcdf.open(url,'NOWRITE')
dimid = netcdf.inqDimID(ncid,'NST')
netcdf.inqDim(ncid,dimid)
[dimname,NST]=netcdf.inqDim(ncid,dimid)
dimid = netcdf.inqDimID(ncid,'Nbed')
netcdf.inqDim(ncid,dimid)
[dimname,Nbed]=netcdf.inqDim(ncid,dimid)
netcdf.close(ncid)
% some cases have extra non-depositing non-cohesive sediment and a few
% classes that interact with the bed
if cas == 66
NNN = 15 % these don't settle
NND = 4 % these interact with the bed
NST = NST-(NNN+NND)
end
%% Read in basic info
ocean_time = ncread(url,'ocean_time');
h = ncread(url,'h');
%lonr = ncread(url,'lon_rho');
%latr = ncread(url,'lat_rho');
x_rho = ncread(url,'x_rho');whos
y_rho = ncread(url,'y_rho');
pn = ncread(url,'pn');
pm = ncread(url,'pm');
%ang = ncread(url,'angle');
Vtransform = ncread(url,'Vtransform')
Vstretching = ncread(url,'Vstretching')
theta_s = ncread(url,'theta_s')
theta_b = ncread(url,'theta_b')
Tcline = ncread(url,'Tcline')
hc = ncread(url,'hc')
s_rho = ncread(url,'s_rho')
s_w = ncread(url,'s_w')
nt = length(ocean_time);
nz = length(s_rho);
%% Read in sediment sizes
fdiam = ncread(url,'Sd50')
ws = ncread(url,'Wsed')
rhos = ncread(url,'Srho')
tau_ce = ncread(url,'tau_ce')
tau_cd = ncread(url,'tau_cd')
%% this is a 1D run, so use i=3 and j =4
i=3; j=4;
h = squeeze(ncread(url,'h',[i j],[1 1]));
dzr = diff(s_w)
ubar = squeeze(ncread(url,'ubar',[i j 1],[1 1 Inf]));
%zeta = squeeze(ncread(url,'zeta'));
zeta = squeeze(ncread(url,'zeta',[i j 1],[1 1 Inf]));
% vbar is tiny...ignore
% vbar = squeeze(ncread(url,'vbar',[i j 1],[1 1 Inf]));
dzr = diff(s_w);
hw = zeta-h;
%% get time series of depths at one location
%igrid = 1 % rho points
%igrid = 3 % upoints
%igrid = 4 % vpoints
%igrid = 5 % w-velocity points
z_r = NaN*ones(nz,nt);
z_w = NaN*ones(nz+1,nz);
for n=1:nt
z_r(:,n)=squeeze(set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, nz, ...
                1, h, zeta(n),0))';
z_w(:,n)=squeeze(set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, nz, ...
                 5, h, zeta(n),0))';
end
dzw = diff(z_w,1,1);
tz = repmat(ocean_time',50,1);
%% read 15 mud classes, one cell only, for 1D results
nt = length(ocean_time);
nz = length(s_rho);

m = zeros( NST, nz, nt);
for n=1:NST
   ncname = sprintf('mud_%02d',n)
   m(n,:,:)=squeeze(ncread(url,ncname,[i j 1 1],[1 1 Inf Inf]));
end
muds = squeeze(sum(m));
summuds = sum( muds.*dzw);
% total conc
s2d = 1. /(3600.*24.)
figure(1); clf
pcolorjw( s2d*tz, h+z_w, log10(muds+eps))
caxis([-2 2])
colorbar
title('log_{10} Total Floc Concentration')

% fraction-weighted ws
wst = squeeze(sum(repmat(ws(1:NST),1,nz,nt).*m)./sum(m));
figure(2); clf
pcolorjw( s2d*tz, h+z_w, 1e3*wst)
colorbar
title('Floc Settling Velocity (mm/s)')

% fraction-weighted size
figure(3); clf
fdiamt = squeeze(sum(repmat(fdiam(1:NST),1,nz,nt).*m)./sum(m));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamt)
colorbar
title('Floc Diameter (\mum)')

% acoustic response 
Dfv = fdiam(1:NST);
rhofv = rhos(1:NST);
v = zeros(nz,nt);
for jj=1:nt
   for ii=1:nz
      mv = squeeze(m(:,ii,jj));
      v(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,3e6);
   end
end
figure(4)
pcolorjw( s2d*tz, h+z_w, v)
colorbar
title('Acoustic Response')
%% read 15 non-depositing classes, one cell only, for 1D results
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
wst = squeeze(sum(repmat(ws(NST+1:NST+NNN),1,nz,nt).*snn)./sum(snn));
figure(6); clf
pcolorjw( s2d*tz, h+z_w, 1e3*wst)
colorbar
title('Non-depositing Sand Settling Velocity (mm/s)')

% fraction-weighted size
figure(7); clf
fdiamt = squeeze(sum(repmat(fdiam(NST+1:NST+NNN),1,nz,nt).*snn)./sum(snn));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamt)
colorbar
title('Non-depositing Sand Diameter (\mum)')

% acoustic response 
Dfv = fdiam(NST+1:NST+NNN);
rhofv = rhos(NST+1:NST+NNN);
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
%% read 4 depositing classes, one cell only, for 1D results
snd = zeros( NND, nz, nt);
ic = 0;
for n=NNN+1:NNN+NND
   ic = ic+1
   ncname = sprintf('sand_%02d',n)
   snd(ic,:,:)=squeeze(ncread(url,ncname,[i j 1 1],[1 1 Inf Inf]));
end
snds = squeeze(sum(snd));
sumsnds = sum( snds.*dzw);

% total conc
s2d = 1. /(3600.*24.)
figure(9); clf
pcolorjw( s2d*tz, h+z_w, log10(snds+eps))
caxis([-2 2])
colorbar
title('Sand log_{10} Total Concentration')

% fraction-weighted ws
wssdt = squeeze(sum(repmat(ws(NNN+NST+1:NST+NNN+NND),1,nz,nt).*snd)./sum(snd));
figure(10); clf
pcolorjw( s2d*tz, h+z_w, 1e3*wssdt)
colorbar
title('Sand Settling Velocity (mm/s)')

% fraction-weighted size
figure(11); clf
fdiamsdt = squeeze(sum(repmat(fdiam(NNN+NST+1:NST+NNN+NND),1,nz,nt).*snd)./sum(snd));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamsdt)
colorbar
title('Sand Diameter (\mum)')

% acoustic response 
Dfv = fdiam(NNN+NST+1:NST+NNN+NND);
rhofv = rhos(NNN+NST+1:NST+NNN+NND);
vsn = zeros(nz,nt);
for jj=1:nt
   for ii=1:nz
      mv = squeeze(snd(:,ii,jj));
      vsn(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,3e6);
   end
end
figure(12)
pcolorjw( s2d*tz, h+z_w, vsn)
colorbar
title('Sand Acoustic Response')

%% floc + sand combined acoustic response 
Dfv = [fdiam(1:NST); fdiam(NNN+NST+1:NST+NNN+NND) ];
rhofv = [rhos(1:NST); rhos(NNN+NST+1:NST+NNN+NND) ];
vcomb = zeros(nz,nt);
for jj=1:nt
   for ii=1:nz
      mv = [squeeze(m(:,ii,jj)); squeeze(snd(:,ii,jj))];
      vcomb(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,3e6);
   end
end
figure(13)
pcolorjw( s2d*tz, h+z_w, vcomb)
colorbar
title('Combined Acoustic Response')
