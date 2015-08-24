% floc_proc.m - Script to read and plot ROMS .his files


url = 'http://geoport.whoi.edu/thredds/dodsC/peach/data2/aretxabaleta/MVCO/floc/nofloc_test35/ocean_his.nc'
url = 'ocean_his.nc'
ncid = netcdf.open(url,'NOWRITE')
dimid = netcdf.inqDimID(ncid,'NST')
netcdf.inqDim(ncid,dimid)
[dimname,NST]=netcdf.inqDim(ncid,dimid)
dimid = netcdf.inqDimID(ncid,'NBED')
netcdf.inqDim(ncid,dimid)
[dimname,NBED]=netcdf.inqDim(ncid,dimid)
netcdf.close(ncid)
%%
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
NST = 15;
NBED = 3;
%% this is a 1D run, so use i=3 and j =4
i=3; j=4;
h = squeeze(ncread(url,'h',[i j],[Inf Inf]));
dzr = diff(s_w)
ubar = squeeze(ncread(url,'ubar',[i j 1],[1 1 Inf]));
%zeta = squeeze(ncread(url,'zeta'));
zeta = squeeze(ncread(url,'zeta',[i j 1],[1 1 Inf]));
% vbar is tiny...ignore
% vbar = squeeze(ncread(url,'vbar',[i j 1],[1 1 Inf]));
dzr = diff(s_w);
hw = zeta-h;
%% read 15 mud classes and mudmass in the bed
% TODO: complete this to make an airtight check on mass conservation
if(0)

m = zeros( nsed, 8, 7, nz, nt);
bm = zeros( nsed, 8, 7, nb, nt);
mst = zeros(8,7,nt)
for n=1:nsed
   ncname = sprintf('mud_%02d',n)
   ncname2 = sprintf('mudmass_%02d',n)
   %m(n,:,:)=squeeze(ncread(url,ncname,[i j 1 1],[1 1 Inf Inf]));
   m(n,:,:,:,:)=squeeze(ncread(url,ncname));
   bm(n,:,:,:,:)=squeeze(ncread(url,ncname2));
end
end
%% read 15 mud classes, one cell only, for 1D results
nt = length(ocean_time);
nz = length(s_rho);
NST = 15; 
NBED = 3;
m = zeros( NST, nz, nt);
for n=1:NST
   ncname = sprintf('mud_%02d',n)
   m(n,:,:)=squeeze(ncread(url,ncname,[i j 1 1],[1 1 Inf Inf]));
end