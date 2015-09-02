% check_stress - script to compare modeled v. observed stresses
% floc_proc.m - Script to read and plot ROMS .his files
clear

cas = 66
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
tdays = ocean_time/(3600*24);
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
%% other basic stuff
i=3; j=4;
%% compare with observations
rho0 = ncread(url,'rho0')
ustrc = (squeeze(ncread(url,'bustrc',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrc',[i j 1],[1 1 Inf])).^2).^(.5);
ustrcwm = (squeeze(ncread(url,'bustrcwmax',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrcwmax',[i j 1],[1 1 Inf])).^2).^(.5);
load ustar_av
%% time series comparison
figure(1); clf
h1=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.ustrc.^2);
set(h1,'linewidth',3,'color',[.3 .3 .3]);
hold on
h2=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.us.^2);
set(h2,'linewidth',3,'color',[.5 .2 .2]);
h3=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.ustrr.^2);
set(h3,'linewidth',3,'color',[.3 .3 .3]);
h4=plot(tdays,ustrc);
set(h4,'linewidth',3,'color',[.5 .5 .8]);
h5=plot(tdays,ustrcwm);
set(h5,'linewidth',3,'color',[.5 .5 .8]);
ylabel('Stress (Pa)')
%% target diagram
figure(2); clf
[RMSD_star,BIAS,Rustrc]=target_diagram( ustrc, rho0*ustar_av.ustrc.^2, 1, [.2 .2 .9] );
%figure(3); clf
[RMSD_star,BIAS,Rustrcw]=target_diagram( ustrcwm, rho0*ustar_av.ustrr.^2, 1 ,[.9 .2 .2]);
figure(4); clf
plot([0 2],[0 2],'--k')
hold on
h1=plot(rho0*ustar_av.ustrr.^2,ustrcwm,'o')
set(h1,'markersize',10,'color',[.2 .2 .2],'markerfacecolor',[.3 .3 .7])
h2=plot(rho0*ustar_av.ustrc.^2,ustrc,'o')
set(h2,'markersize',10,'color',[.2 .2 .2],'markerfacecolor',[.7 .3 .3])
axis([0. 2. 0 2])
axis 'square'
lt2 = ['\it{\tau_{*c} r^2}=',sprintf('%6.3f',Rustrc)]
lt1 = ['\it{\tau_{*cw} r^2}=',sprintf('%6.3f',Rustrcw)]
legend([h2; h1],lt2,lt1)
ts = sprintf('Run %d',cas)
xlabel('Data')
ylabel('Model')

% check lags