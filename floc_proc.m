% floc_proc.m - Script to read and plot ROMS .his files
clear
cas = 78
fn = 0;
iplot = 0; % if true, save plots
%url = sprintf('ocean_his%2d.nc',cas)
% url = 'http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/aretxabaleta/MVCO/ocean_his_44.nc'
url = sprintf('http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/aretxabaleta/MVCO/ocean_his_%02d.nc', cas)
% Read in NST and Nbed instead of loading in huge files to infer their size
ncid = netcdf.open(url,'NOWRITE');
dimid = netcdf.inqDimID(ncid,'NST')
netcdf.inqDim(ncid,dimid);
[dimname,NST]=netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'Nbed');
netcdf.inqDim(ncid,dimid);
[dimname,Nbed]=netcdf.inqDim(ncid,dimid)
netcdf.close(ncid);
% some cases have extra non-depositing non-cohesive sediment and a few
% classes that interact with the bed
if cas >= 66
    NNN = 15 % these don't interact with the bed
    NND = 4  % these interact with the bed
    NCS = NST-(NNN+NND)
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
s_rho = ncread(url,'s_rho');
s_w = ncread(url,'s_w');
nt = length(ocean_time)
nz = length(s_rho)
%% Read in sediment sizes
fdiam = ncread(url,'Sd50');
ws = ncread(url,'Wsed');
rhos = ncread(url,'Srho');
tau_ce = ncread(url,'tau_ce');
tau_cd = ncread(url,'tau_cd');
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
elev = h+z_w(:,1);
%% read 15 mud classes, one cell only, for 1D results
nt = length(ocean_time);
nz = length(s_rho);

m = zeros( NCS, nz, nt);
for n=1:NCS
    ncname = sprintf('mud_%02d',n)
    m(n,:,:)=squeeze(ncread(url,ncname,[i j 1 1],[1 1 Inf Inf]));
end
muds = squeeze(sum(m));
summuds = sum( muds.*dzw);
%% acoustic response of flocs - v1, v25, v4
Dfv = fdiam(1:NCS);
rhofv = rhos(1:NCS);
v1 = zeros(nz,nt);
v25 = zeros(nz,nt);
v4 = zeros(nz,nt);
for jj=1:nt
    for ii=1:nz
        mv = squeeze(m(:,ii,jj));
        v1(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,1e6);
        v25(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,2.5e6);
        v4(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,4e6);
    end
end
%% optical response to flocs of ac9 and LISST - r_ac9f, r_lisstf
r_ac9f=zeros(nz,nt);
r_lisstf=zeros(nz,nt);
if(exist('isfloc','var')~=1),isfloc=1; end
% import D, nf, c_ac9, c_mass, c_lisst, and same for sand
if(exist('c_ac9','var')~=1),
    load c_coeffs_crs
end
% skip next step because nf=2 is column 1
% c_ac9_nf = interp1(nf_optics',c_ac9',fnf)';
c_ac9_i = interp1(D_optics*1e-6',c_ac9(:,1),Dfv')';
c_lisst_i = interp1(D_optics*1e-6',c_lisst(:,1),Dfv')';
r_ac9=sum(mv.*c_ac9_i);
for jj=1:nt
    for ii=1:nz
        mv = squeeze(m(:,ii,jj));
        r_ac9f(ii,jj)=sum(mv.*c_ac9_i);
        r_lisstf(ii,jj)=sum(mv.*c_lisst_i);
    end
end
%% read 15 non-depositing classes, one cell only, for 1D results
if(0)
    load_non
end
%% read 4 depositing sand classes, one cell only, for 1D results
snd = zeros( NND, nz, nt);
ic = 0;
for n=NNN+1:NNN+NND
    ic = ic+1
    ncname = sprintf('sand_%02d',n)
    snd(ic,:,:)=squeeze(ncread(url,ncname,[i j 1 1],[1 1 Inf Inf]));
end
snds = squeeze(sum(snd));
sumsnds = sum( snds.*dzw);
%% acoustic response of sand - vsn1, vsn25, and vsn4 (three frequencies)
Dfv = fdiam(NNN+NCS+1:NCS+NNN+NND);
rhofv = rhos(NNN+NCS+1:NCS+NNN+NND);
vsn1 = zeros(nz,nt);
vsn25 = zeros(nz,nt);
vsn3 = zeros(nz,nt);
for jj=1:nt
    for ii=1:nz
        mv = squeeze(snd(:,ii,jj));
        vsn1(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,1e6);
        vsn25(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,2.5e6);
        vsn4(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,4e6);
    end
end
%% optical response of sand to ac9 and LISST - r_ac9s, r_lissts
Dfv = fdiam(NNN+NCS+1:NCS+NNN+NND);
r_ac9s=zeros(nz,nt);
r_lissts=zeros(nz,nt);
% import D, nf, c_ac9, c_mass, c_lisst, and same for sand
if(exist('c_ac9_sand','var')~=1),
    load c_coeffs_crs
end
% skip next step because nf=2 is column 1
% c_ac9_nf = interp1(nf_optics',c_ac9',fnf)';
c_ac9_is = interp1(D_optics*1e-6',c_ac9_sand(:),Dfv')';
c_lisst_is = interp1(D_optics*1e-6',c_lisst_sand(:),Dfv')';
for jj=1:nt
    for ii=1:nz
        mv = squeeze(snd(:,ii,jj));
        r_ac9s(ii,jj)=sum(mv.*c_ac9_is);
        r_lissts(ii,jj)=sum(mv.*c_ac9_is);
    end
end
%% floc + sand combined acoustic response - vcomb
Dfv = [fdiam(1:NCS); fdiam(NNN+NCS+1:NCS+NNN+NND) ];
rhofv = [rhos(1:NCS); rhos(NNN+NCS+1:NCS+NNN+NND) ];
vcomb1 = zeros(nz,nt);
vcomb25 = zeros(nz,nt);
vcomb4 = zeros(nz,nt);
for jj=1:nt
    for ii=1:nz
        mv = [squeeze(m(:,ii,jj)); squeeze(snd(:,ii,jj))];
        vcomb1(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,1e6);
        vcomb25(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,2.5e6);
        vcomb4(ii,jj) = acoustics(Dfv(:),mv(:),rhofv(:),.2,4e6);
    end
end
%% fraction-weighted size of both sand and flocs - fdiamall
fdiamall = squeeze((sum(repmat(fdiam(1:NCS),1,nz,nt).*m)...
    +        sum(repmat(fdiam(NNN+NCS+1:NCS+NNN+NND),1,nz,nt).*snd))...
    ./(sum(snd)+sum(m)));

