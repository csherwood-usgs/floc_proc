% floc_proc.m - Script to read and plot ROMS .his files
clear

cas = 50
%url = sprintf('ocean_his%2d.nc',cas)
% url = 'http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/aretxabaleta/MVCO/ocean_his_44.nc'
url = sprintf('http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/aretxabaleta/MVCO/ocean_his_%02d.nc', cas)
% Read in NST and Nbed instead of loading in huge files to infer their size
ncid = netcdf.open(url,'NOWRITE')
dimid = netcdf.inqDimID(ncid,'NST')
netcdf.inqDim(ncid,dimid)
[dimname,NST]=netcdf.inqDim(ncid,dimid)
dimid = netcdf.inqDimID(ncid,'Nbed')
netcdf.inqDim(ncid,dimid)
[dimname,Nbed]=netcdf.inqDim(ncid,dimid)
netcdf.close(ncid)
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
%% read 15 mud classes and mudmass in the bed
% TODO: complete this to make an airtight check on mass conservation
% (use evaluate_sediment.m to do this)

%% read 15 mud classes, one cell only, for 1D results
nt = length(ocean_time);
nz = length(s_rho);
NST = 15; 
Nbed = 3;
m = zeros( NST, nz, nt);
for n=1:NST
   ncname = sprintf('mud_%02d',n)
   m(n,:,:)=squeeze(ncread(url,ncname,[i j 1 1],[1 1 Inf Inf]));
end
muds = squeeze(sum(m));
summuds = sum( muds.*dzw);
%% total conc
s2d = 1. /(3600.*24.)
figure(1); clf
pcolorjw( s2d*tz, h+z_w, log10(muds+eps))
colorbar
title('log_{10} Total Concentration')
% fraction-weighted ws
wst = squeeze(sum(repmat(ws,1,nz,nt).*m)./sum(m));
figure(2); clf
pcolorjw( s2d*tz, h+z_w, 1e3*wst)
colorbar
title('Settling Velocity (mm/s)')
% fraction-weighted size
figure(3); clf
fdiamt = squeeze(sum(repmat(fdiam,1,nz,nt).*m)./sum(m));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamt)
colorbar
title('Diameter (\mum)')
%% compare with observations
rho0 = ncread(url,'rho0')
ustrc = (squeeze(ncread(url,'bustrc',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrc',[i j 1],[1 1 Inf])).^2).^(.5);
ustrcwm = (squeeze(ncread(url,'bustrcwmax',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrcwmax',[i j 1],[1 1 Inf])).^2).^(.5);
load ustar_av
%%
figure(4); clf
h1=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.ustrc.^2);
set(h1,'linewidth',3,'color',[.3 .3 .3]);
hold on
h2=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.us.^2);
set(h2,'linewidth',3,'color',[.2 .2 .2]);
h3=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.ustrr.^2);
set(h3,'linewidth',3,'color',[.3 .3 .3]);
h4=plot(s2d*tz,ustrc);
set(h4,'linewidth',3,'color',[.4 .4 .6]);
h5=plot(s2d*tz,ustrcwm);
set(h5,'linewidth',3,'color',[.4 .4 .6]);
ylabel('Stress (Pa)')