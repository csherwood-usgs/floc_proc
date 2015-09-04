% floc_proc.m - Script to read and plot ROMS .his files
clear
cas = 70
iplot = 0; % if true, save plots
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
if cas >= 66
   NNN = 15 % these don't settle
   NND = 4  % these interact with the bed
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

%% floc conc plots
s2d = 1. /(3600.*24.)
figure(1); clf
subplot(311)
pcolorjw( s2d*tz, h+z_w, log10(muds+eps))
caxis([-2.5 1.5])
set(gca,'xticklabel','','fontsize',14)
colorbar
title('log_{10} Total Floc Concentration (kg/m^3)')

% fraction-weighted ws
wst = squeeze(sum(repmat(ws(1:NST),1,nz,nt).*m)./sum(m));
subplot(312)
pcolorjw( s2d*tz, h+z_w, 1e3*wst)
set(gca,'xticklabel','','fontsize',14)
colorbar
title('Floc Settling Velocity (mm/s)')
ylabel('Elevation (m)')

% fraction-weighted size
subplot(313)
fdiamt = squeeze(sum(repmat(fdiam(1:NST),1,nz,nt).*m)./sum(m));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamt)
colorbar
title('Floc Diameter (\mum)')
set(gca,'fontsize',14)
xlabel('Days','fontsize',16)
pfn=sprintf('floc_conc_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% acoustic response
Dfv = fdiam(1:NST);
rhofv = rhos(1:NST);
v = zeros(nz,nt);
for jj=1:nt
   for ii=1:nz
      mv = squeeze(m(:,ii,jj));
      v(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,3e6);
   end
end
%% flocs optical response
r_ac9f=zeros(nz,nt);
if(exist('isfloc','var')~=1),isfloc=1; end
% import D, nf, c_ac9, c_mass, c_lisst, and same for sand
if(exist('c_ac9','var')~=1),
   load c_coeffs_crs
end
% skip next step because nf=2 is column 1
% c_ac9_nf = interp1(nf_optics',c_ac9',fnf)';
c_ac9_i = interp1(D_optics*1e-6',c_ac9(:,1),Dfv')';
r_ac9=sum(mv.*c_ac9_i);
for jj=1:nt
   for ii=1:nz
      mv = squeeze(m(:,ii,jj));
      r_ac9f(ii,jj)=sum(mv.*c_ac9_i);
   end
end
%% read 15 non-depositing classes, one cell only, for 1D results
if(0)
   
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
   xlabel('Days')
   ylabel('Elevation (m)')
   xlim([0 32.8])
   caxis([ 0 .3])
   ylim([0 3])
end
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
%% sand conc plots
s2d = 1. /(3600.*24.)
figure(2); clf
subplot(311)
pcolorjw( s2d*tz, h+z_w, log10(snds+eps))
set(gca,'xticklabel','','fontsize',14)
caxis([-2.5 1.5])
colorbar
title('log_{10} Total Sand Concentration (kg/m^3)')

% fraction-weighted ws
wssdt = squeeze(sum(repmat(ws(NNN+NST+1:NST+NNN+NND),1,nz,nt).*snd)./sum(snd));
subplot(312)
pcolorjw( s2d*tz, h+z_w, 1e3*wssdt)
set(gca,'xticklabel','','fontsize',14)
colorbar
ylabel('Elevation (m)','fontsize',16)
title('Sand Settling Velocity (mm/s)')

% fraction-weighted size
subplot(313)
fdiamsdt = squeeze(sum(repmat(fdiam(NNN+NST+1:NST+NNN+NND),1,nz,nt).*snd)./sum(snd));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamsdt)
colorbar
set(gca,'fontsize',14)
xlabel('Days','fontsize',16)
title('Sand Diameter (\mum)')
pfn=sprintf('sand_conc_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% acoustic response
Dfv = fdiam(NNN+NST+1:NST+NNN+NND);
rhofv = rhos(NNN+NST+1:NST+NNN+NND);
vsn = zeros(nz,nt);
for jj=1:nt
   for ii=1:nz
      mv = squeeze(snd(:,ii,jj));
      vsn(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,3e6);
   end
end
%% optical response
Dfv = fdiam(NNN+NST+1:NST+NNN+NND);
r_ac9s=zeros(nz,nt);
% import D, nf, c_ac9, c_mass, c_lisst, and same for sand
if(exist('c_ac9_sand','var')~=1),
   load c_coeffs_crs
end
% skip next step because nf=2 is column 1
% c_ac9_nf = interp1(nf_optics',c_ac9',fnf)';
c_ac9_is = interp1(D_optics*1e-6',c_ac9_sand(:),Dfv')';
for jj=1:nt
   for ii=1:nz
      mv = squeeze(snd(:,ii,jj));
      r_ac9s(ii,jj)=sum(mv.*c_ac9_is);
   end
end
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
%% Plot acoustic response
figure(4); clf
subplot(311)
pcolorjw( s2d*tz, h+z_w, v)
colorbar
title('Flocs Acoustic Response')
set(gca,'xticklabel','','fontsize',14)
%xlabel('Days')
%ylabel('Elevation (m)')
xlim([0 32.8])
caxis([ 0 .5])
ylim([0 3])
subplot(312)
pcolorjw( s2d*tz, h+z_w, vsn)
colorbar
title('Sand Acoustic Response')
set(gca,'xticklabel','','fontsize',14)
%xlabel('Days')
ylabel('Elevation (m)','fontsize',16)
xlim([0 32.8])
caxis([ 0 .5])
ylim([0 3])

subplot(313)
pcolorjw( s2d*tz, h+z_w, vcomb)
colorbar
title('Combined Acoustic Response')
set(gca,'fontsize',14)
xlabel('Days','fontsize',16)
%ylabel('Elevation (m)')
xlim([0 32.8])
caxis([ 0 .5])
ylim([0 3])

pfn=sprintf('both_acoustic_response_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% Plot optical responses
figure(5); clf
subplot(311)
pcolorjw( s2d*tz, h+z_w, log10(r_ac9f+eps))
colorbar
title('Flocs Optical (ac-9) Response')
set(gca,'xticklabel','','fontsize',14)
%ylabel('Elevation (m)')
xlim([0 32.8])
caxis([-3 1])
ylim([0 3])

subplot(312)
pcolorjw( s2d*tz, h+z_w, log10(r_ac9s+eps))
colorbar
title('Sand Optical Response')
set(gca,'xticklabel','','fontsize',14)
ylabel('Elevation (m)','fontsize',16)
xlim([0 32.8])
caxis([ -3 1])
ylim([0 3])

subplot(313)
pcolorjw( s2d*tz, h+z_w, log10(r_ac9f+r_ac9s+eps))
colorbar
title('Combined Optical Response')
set(gca,'fontsize',14)
xlabel('Days','fontsize',16)
%ylabel('Elevation (m)','fontsize',14)
xlim([0 32.8])
caxis([ -3 1])
ylim([0 3])

pfn=sprintf('both_optical_response_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%%
load cmap_plusminus
figure(14)
pcolorjw( s2d*tz, h+z_w, vcomb-vsn)
colorbar
title('Difference: Combined - Sand Acoustic Response','fontsize',16)
xlabel('Days','fontsize',16)
ylabel('Elevation (m)','fontsize',16)
xlim([0 32.8])
caxis([ -.25 .25])
ylim([0 3])
pfn=sprintf('acoustic_diff_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% overall acoustic response to conc
tsed = muds+snds;
figure(15); clf
h1=plot(tsed(:),vcomb(:),'.','color',[.3 .3 .3]);
set(h1,'markersize',14)
hold on
h2=plot(muds(:),v(:),'.','color',[.4 0 0]);
set(h2,'markersize',14)
h3=plot(snds(:),vsn(:),'.','color',[1. .4 0]);
set(h3,'markersize',14)
xlabel('Mass Concentration (kg/m^3)','fontsize',16)
ylabel('Acoustic Response','fontsize',16)
h4=legend([h3;h2;h1],'Sand','Flocs','Combined')
set(h4,'fontsize',14)
pfn = sprintf('acoustic_response_scatter%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end