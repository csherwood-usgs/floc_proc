% check_stress - script to compare modeled v. observed stresses
% floc_proc.m - Script to read and plot ROMS .his files
clear

cas = 70
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
tauc = (squeeze(ncread(url,'bustrc',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrc',[i j 1],[1 1 Inf])).^2).^(.5);
taucwm = (squeeze(ncread(url,'bustrcwmax',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrcwmax',[i j 1],[1 1 Inf])).^2).^(.5);
load ustar_av
%% time series comparison
figure(1); clf
h1=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.ustrc.^2);
set(h1,'linewidth',3,'color',[.3 .3 .3]);
hold on
%h2=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.us.^2);
%set(h2,'linewidth',3,'color',[.5 .2 .2]);
h3=plot(ustar_av.dn-ustar_av.dn(1),rho0*ustar_av.ustrr.^2);
set(h3,'linewidth',3,'color',[.3 .3 .3]);
h4=plot(tdays,tauc);
set(h4,'linewidth',3,'color',[.7 .3 .3]);
h5=plot(tdays,taucwm);
set(h5,'linewidth',3,'color',[.3 .3 .7]);
set(gca,'fontsize',14)
ylabel('Stress (Pa)','fontsize',16)
xlabel('Days','fontsize',16)
xlim([0 32.8])
ylim([0 4])
h6=legend([h3;h4;h5],'Obs','$\it{\tau_{*c}}$','$\it{\tau_{*cw}}$')
set(h6,'interpreter','latex','fontsize',16)
%% target diagram
ou=nanmean(abs(rho0*ustar_av.us.^2-rho0*ustar_av.ustrc.^2))/nanstd(rho0*ustar_av.us.^2)
figure(2); clf
[RMSD_star,BIAS,Rustrc]=target_diagram( tauc, rho0*ustar_av.ustrc.^2, 1, [.2 .2 .9],ou );
set(gca,'fontsize',16)
%figure(3); clf
[RMSD_star,BIAS,Rustrcw]=target_diagram( taucwm, rho0*ustar_av.ustrr.^2, 1 ,[.9 .2 .2]);
set(gca,'fontsize',16)
figure(4); clf
plot([0 2],[0 2],'--k')
hold on
h1=plot(rho0*ustar_av.ustrr.^2,taucwm,'o')
set(h1,'markersize',10,'color',[.2 .2 .2],'markerfacecolor',[.3 .3 .7])
h2=plot(rho0*ustar_av.ustrc.^2,tauc,'o')
set(h2,'markersize',10,'color',[.2 .2 .2],'markerfacecolor',[.7 .3 .3])
set(gca,'fontsize',14,'xtick',[0:.5:2],'ytick',[0:.5:2])
axis([0. 2. 0 2])
axis 'square'
lt2 = ['\it{\tau_{*c} r^2}=',sprintf('%6.3f',Rustrc)]
lt1 = ['\it{\tau_{*cw} r^2}=',sprintf('%6.3f',Rustrcw)]
h6=legend([h2; h1],lt2,lt1)
set(h6,'fontsize',16)
ts = sprintf('Run %d',cas)
xlabel('Data')
ylabel('Model')