% floc_proc.m - Script to read and plot ROMS .his files

% disable next three variables if running many cases:
% cas = 107
% fid = 1;
% iplot = 0; % if true, save plots

% convert seconds to days
s2d = 1. /(3600.*24.);

% File with case info
fprintf(fid2,'%d, ',cas)

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
NNN = 0 % these don't interact with the bed
NND = 4  % these interact with the bed
NCS = NST-(NNN+NND)
if cas >= 66 & cas <=93
   NNN = 15 % these don't interact with the bed
   NND = 4  % these interact with the bed
   NCS = NST-(NNN+NND)
end
fprintf(fid2,'%d, %d, ',NCS,NST);
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
fprintf(fid2,'%d, %d, ',nt,nz);
%% Read in sediment sizes
fdiam = ncread(url,'Sd50');
ws = ncread(url,'Wsed');
rhos = ncread(url,'Srho');
tau_ce = ncread(url,'tau_ce');
tau_cd = ncread(url,'tau_cd');
%% this is a 1D run, so use i=3 and j =4
i=3; j=4;
h = squeeze(ncread(url,'h',[i j],[1 1]));
dzr = diff(s_w);
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
tz = repmat(ocean_time',nz,1);
elev = h+z_w(:,1);
mean_dz = mean(dzw(:));
% index for cell closest to 0.5 m
iz50 = find(elev>=(0.5-(0.5*mean_dz)),1,'first');
fprintf(fid2,'%4.2f, %d, %4.2f, ',mean_dz,iz50,elev(iz50));
%% Read in stresses
rho0 = ncread(url,'rho0');
bustrc=squeeze(ncread(url,'bustrc',[i j 1],[1 1 Inf]));
tauc = (bustrc.^2+squeeze(ncread(url,'bvstrc',[i j 1],[1 1 Inf])).^2).^(.5);
taucwm = (squeeze(ncread(url,'bustrcwmax',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrcwmax',[i j 1],[1 1 Inf])).^2).^(.5);
ustrc = sqrt(tauc./rho0);
ustrcw = sqrt(taucwm./rho0);
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
% characterize mud load at 0.5 mab (skip first six hours...startup)
mean_mud_conc = mean(muds(iz50,6:end));
max_mud_conc = max(muds(iz50,6:end));
fprintf(fid2,'%9.6f, %9.6f, ',mean_mud_conc,max_mud_conc);
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
% characterize sand load (skip first six hours...startup)
mean_sand_conc = mean(snds(iz50,6:end));
max_sand_conc = max(snds(iz50,6:end));
fprintf(fid2,'%9.6f, %9.6f, ',mean_sand_conc,max_sand_conc);
%% v1, v25, v4 - acoustic response to flocs
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
%% vsn1, vsn25, and vsn4 - acoustic response to sand
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
%% vcomb1, vcomb25, vcomb4 - floc + sand combined acoustic response
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
%%  r_ac9f, r_lisstf - optical response to flocs of ac9 and LISST
Dfv = fdiam(1:NCS);
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
% zero out LISST sizes that are too big
c_lisst_i(Dfv>500.)=0;
for jj=1:nt
   for ii=1:nz
      mv = squeeze(m(:,ii,jj));
      r_ac9f(ii,jj)=sum(mv.*c_ac9_i);
      r_lisstf(ii,jj)=sum(mv.*c_lisst_i);
   end
end
%% r_ac9s, r_lissts - optical response to sand of ac9 and LISST
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
% zero out LISST sizes that are too big
c_lisst_is(Dfv>500.)=0;
for jj=1:nt
   for ii=1:nz
      mv = squeeze(snd(:,ii,jj));
      r_ac9s(ii,jj)=sum(mv.*c_ac9_is);
      r_lissts(ii,jj)=sum(mv.*c_ac9_is);
   end
end
%% r_ac9c, r_lisstc - optical response of LISST  to sand and flocs combined
r_ac9c=zeros(nz,nt);
r_lisstc=zeros(nz,nt);
for jj=1:nt
   for ii=1:nz
      mv = [c_ac9_i.*squeeze(m(:,ii,jj)); c_ac9_is.*squeeze(snd(:,ii,jj))];
      r_ac9c(ii,jj)=sum(mv);
      mv = [c_lisst_i.*squeeze(m(:,ii,jj)); c_lisst_is.*squeeze(snd(:,ii,jj))];
      r_lisstc(ii,jj)=sum(mv);
   end
end
stats(r_ac9c(3,:));
stats(r_lisstc(3,:));
%% calculate the ac9 posteriori calibration coefficient for the ac9 at 1 mab
% this is at circa one meter...and gets overwritten for all depths in the
% next section
model_mass_conc = squeeze( sum(m(:,iz50,:)) + sum(snd(:,iz50,:)));
ac9_response = r_ac9c(iz50,:);
abss25_response = vcomb25(iz50,:);

model_color = [.8 .2 .2];
data_color = [.3 .3 .3];
figure(1); clf;
% ac9
subplot(221)
[ac9_offset,ac9_slope,ac9_r2, sa, sb, hdot ] = lsfit(model_mass_conc(model_mass_conc<.4),ac9_response(model_mass_conc<.4)',1)
xlim([0. .25])
title('ac9')
xlabel('Modeled Concentration (kg/m^3)')
ylabel('Modeled ac9 Response')
ts = sprintf('%8.4f x C + %8.4f\nr^2 = %4.2f',ac9_slope,ac9_offset,ac9_r2)
ax = axis;
ht = text(.05*ax(2),.9*ax(4),ts);
fprintf(fid2,'%8.4f, %8.4f, %4.2f',ac9_slope,ac9_offset,ac9_r2)

subplot(223)
h1=plot(s2d*ocean_time, model_mass_conc,'linewidth',2,'color',data_color);
hold on
h2=plot(s2d*ocean_time,(ac9_response-ac9_offset)/ac9_slope,'linewidth',2,'color',model_color );
ylim([0 .25])
ylabel('Concentration (kg/m^3)')
legend([h1;h2],'Modeled Concentration','Modeled Response','location','southeast')
xlabel('Time (days)')

% ABSS
subplot(222)
[abss_offset,abss_slope,abss_r2, sa, sb, hdot ] = lsfit(model_mass_conc(model_mass_conc<.4),abss25_response(model_mass_conc<.4)',1)
xlim([0. .25])
title('2.5 MHz ABSS Response')
xlabel('Modeled Concentration (kg/m^3)')
ylabel('Modeled ABSS Response')
ts = sprintf('%8.4f x C + %8.4f\nr^2 = %4.2f',abss_slope,abss_offset,abss_r2)
ax = axis;
ht = text(.05*ax(2),.9*ax(4),ts);
fprintf(fid2,'%8.4f, %8.4f, %4.2f',abss_slope,abss_offset,abss_r2)

subplot(224)
h1=plot(s2d*ocean_time, model_mass_conc,'linewidth',2,'color',data_color);
hold on
h2=plot(s2d*ocean_time,(abss25_response-abss_offset)/abss_slope,'linewidth',2,'color',model_color );
%h2=plot(s2d*ocean_time,(abss25_response)/abss_slope,'linewidth',2,'color',model_color );
legend([h1;h2],'Modeled Concentration','Modeled Response','location','southeast')
ylim([0 .25])
ylabel('Concentration (kg/m^3)')
xlabel('Time (days)')

pfn = sprintf('calib_plot_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%%
model_mass_conc = squeeze( sum(m(:,:,:)) + sum(snd(:,:,:)));
model_ac9_response = (r_ac9c(:,:)-ac9_offset)/ac9_slope;
% include the offset
model_abss25_response = (vcomb25(:,:)-abss_offset)/abss_slope; 
% ...or neglect the offset
%model_abss25_response = (vcomb25(:,:))/abss_slope;
%% fdiamall - mass fraction-weighted size of both sand and flocs - fdiamall
fdiamall = squeeze((sum(repmat(fdiam(1:NCS),1,nz,nt).*m)...
   +        sum(repmat(fdiam(NNN+NCS+1:NCS+NNN+NND),1,nz,nt).*snd))...
   ./(sum(snd)+sum(m)));
all_d50 = 1e6*median(fdiamall(iz50,6:end))
%% fdiamlisst - volume fraction-weighted size, omit diameters outside LISST range - fdiamlisst
% calc weighting for convertion to volume concentrations
fvwht=1. ./rhos(1:NCS)
svwht= 1. ./rhos(NNN+NCS+1:NCS+NNN+NND)
% copy mass concentrations
ml = m(1:NCS,:,:).*repmat(fvwht,1,nz,nt);
sndl = snd.*repmat(svwht,1,nz,nt);

% elimnate sizes outside LISST range
toobig = find(fdiam(1:NCS)>500e-6);
toosmall = find(fdiam(1:NCS)<2.5e-6);
ml(toobig,:,:)=0;
ml(toosmall,:,:)=0;
toobigs =find(fdiam(NNN+NCS+1:NCS+NNN+NND)>500e-6);
toosmalls =find(fdiam(NNN+NCS+1:NCS+NNN+NND)<2.5e-6);
sndl(toobigs,:,:)=0;
sndl(toosmalls,:,:)=0;
fdiamlisst = squeeze((sum(repmat(fdiam(1:NCS),1,nz,nt).*ml)...
   +                 sum(repmat(fdiam(NNN+NCS+1:NCS+NNN+NND),1,nz,nt).*sndl))) ...
   ./squeeze(sum(sndl)+sum(ml));
lisst_d50 = 1e6*median(fdiamlisst(iz50,6:end))
fprintf(fid2,'%6.1f, %6.1f, ',all_d50,lisst_d50)
%% calculate apparent settling velocity
ic = 1;
izfirst = find(elev>=0.3,1,'first');
izlast = find(elev<=1.9,1,'last');
iz = izfirst:izlast;
za = 0.1;
z=elev(iz);
fprintf(1,'Elevations for pfit between 0.3 and 1.9:\n');
for ii=1:length(iz),fprintf(1,'  %d  %f\n',iz(ii),z(ii));,end

clear pfa pfs pfm pfv pfac9 pfalisst
for ii=1:nt
   figure(2); clf
   t = s2d*ocean_time(ii);
   c = squeeze(muds(iz,ii)+snds(iz,ii));
   hold on
   pf = pfit( c, z, 1, za);
   c = squeeze(snds(iz,ii));
   pfs = pfit( c, z, 1, za);
   c = squeeze(muds(iz,ii));
   pfm = pfit( c, z, 1, za);
   
   c = squeeze(model_abss25_response(iz,ii));
   pfvcombc = pfit( c, z, 0, za);
   c = squeeze(model_ac9_response(iz,ii));
   pfac9c = pfit( c, z, 0, za);
   c = squeeze(r_lisstc(iz,ii));
   pfalisstc = pfit( c, z, 0, za);
   
   ylim( [.1 2.5] );
   xlim( [.01 1] );
   ylabel('Elevation (mab)')
   xlabel('Concentration ( m^{-1} )');
   
   ts = sprintf('%6.2f\nN=%d\nCa=%7.2f\np=% 5.2f\nr^2=%06.4f\n',...
      t,pf.N,pf.Ca,-pf.p,pf.r2);
   text(.2,1.5,ts)
   shg
   %pause(.2)
   pfa(ic)=pf;
   pfsa(ic)=pfs;
   pfma(ic)=pfm;
   pfv(ic)=pfvcombc;
   pfac9(ic)=pfac9c;
   pfalisst(ic)=pfalisstc;
   ic = ic+1;
end
%% plot the pfit results
figure(3); clf
subplot(311)
plot(s2d*ocean_time,zeros(size(ocean_time)),'--k');
hold on
h1=plot(s2d*ocean_time,sign(bustrc).*ustrc,'linewidth',2,'color',[.1 .1 .4]);
h2=plot(s2d*ocean_time,ustrcw,'linewidth',2,'color',[.3 .3 .6]);
h3=plot(ydd-18./24.,ustrcw2h,'linewidth',2,'color',data_color);
ylabel('{\it{u*, u*_{cw}}} (m/s)')
legend([h3;h2;h1],'M94 {\itu*_{cw}}','Model {\itu*_{cw}}','Model {\itu*_c}')

subplot(312)
h1=plot(s2d*ocean_time,-1e3*[pfa.p]'.*(0.41*ustrc),'linewidth',2);
hold on
h2=plot(s2d*ocean_time,-1e3*[pfsa.p]'.*(0.41*ustrc),'linewidth',2);
h3=plot(s2d*ocean_time,-1e3*[pfma.p]'.*(0.41*ustrc),'linewidth',2);
h4=plot(s2d*ocean_time,-1e3*[pfv.p]'.*(0.41*ustrc),'linewidth',2);
h5=plot(s2d*ocean_time,-1e3*[pfac9.p]'.*(0.41*ustrc),'linewidth',2);
h6=plot(s2d*ocean_time,-1e3*[pfalisst.p]'.*(0.41*ustrc),'linewidth',2);
legend([h1;h2;h3;h4;h5;h6],'All','Sand','Flocs','2.5 MHz','ac9','LISST');
ylabel('{\it{w_s}} (mm/s)')

if(0) % raw slopes...need this?
   subplot(413)
h1=plot(s2d*ocean_time,[pfa.p],'linewidth',2);
hold on
h2=plot(s2d*ocean_time,[pfv.p],'linewidth',2);
h3=plot(s2d*ocean_time,[pfac9.p],'linewidth',2);
h4=plot(s2d*ocean_time,[pfalisst.p],'linewidth',2);
legend([h1;h2;h3;h4],'Mass','2.5 MHz','ac9','LISST');
end

subplot(313)
h1=plot(s2d*ocean_time,[pfa.r2],'linewidth',2);
hold on
h2=plot(s2d*ocean_time,[pfv.r2],'linewidth',2);
h3=plot(s2d*ocean_time,[pfac9.r2],'linewidth',2);
h4=plot(s2d*ocean_time,[pfalisst.r2],'linewidth',2);
legend([h1;h2;h3;h4],'Mass','2.5 MHz','ac9','LISST');
pfn = sprintf('pfit_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
ylabel('{\itr}^2')
ylim([.5 1])

pfn = sprintf('ustrcw_ws_r2_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% compare with data - u*, u*cw
% 
% add .7 days to model time
if(0)
figure(4); clf
h6=plot( (18./24.)+s2d*ocean_time,rho0*ustrcw.^2,'linewidth',2);
hold on
h2=plot(ydd,rho0*ustrcw2h.^2,'linewidth',2);
pfn = sprintf('ustrcw_data_model_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
end
%% ABSS ws
figure(5); clf
ok = find(abss2_r22h>0.8);
h1=plot(ydd(ok),1000*abss2_ws2h(ok),'linewidth',2,'color',data_color);
hold on
h2=plot( (18./24.)+s2d*ocean_time,-1e3*[pfv.p]'.*(0.41*ustrc),'linewidth',2,'color',model_color);
title('ABSS Settling Velocity')
ylabel('{\it{w_s}} (mm/s)')
legend([h1;h2],'Data','Model')
pfn = sprintf('abss_ws_data_model_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% LISST ws
figure(6); clf
h2=plot( (18./24.)+s2d*ocean_time,-1e3*[pfalisst.p]'.*(0.41*ustrc),'linewidth',2,'color',model_color);
hold on 
ok = find(LISSTattn_r22h>0.8);
h1=plot(ydd(ok),1000*LISSTattn_ws2h(ok),'linewidth',2,'color',data_color);
title('LISST Settling Velocity')
ylabel('{\it{w_s}} (mm/s)')
legend([h2;h1],'Model','Data')
pfn = sprintf('LISST_ws_data_model_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% LISST D50
figure(7); clf
h2=plot( (18./24.)+s2d*ocean_time,1e6*fdiamlisst(10,:),'linewidth',2,'color',model_color*1.2);
hold on
h3=plot( (18./24.)+s2d*ocean_time,1e6*fdiamlisst(2,:),'linewidth',2,'color',model_color);
%h4=plot( (18./24.)+s2d*ocean_time,1e6*fdiamlisst(30,:),'linewidth',2,'color',model_color);
hold on
h1=plot(ydd,NX(:,9),'linewidth',2,'color',data_color);
title('LISST D50')
ylabel('{\itD}_{50} (\mum)')
xlabel('Time (days)')
legend([h1; h3],'Data','Model')
pfn = sprintf('LISST_D50_data_model_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% ABSS Response
% model elev(2) is 0.23 m, elevl(8) is 1.6 ma
% data bin(1) is .15 m, bin(9) is 1.72
% 15 is an calibration constant...could be, say, the system constant for
% the modeled ABSS values

figure(8); clf
h3=plot(ydd,NX(:,24),'linewidth',2,'color',data_color);
hold on
h4=plot(ydd,NX(:,27),'linewidth',2,'color',data_color*1.5);
set(h4,'color',[.4 .4 .4])
h1=plot( (18./24.)+s2d*ocean_time,model_abss25_response(8,:),'linewidth',2,'color',model_color);
h2=plot( (18./24.)+s2d*ocean_time,(model_abss25_response(2,:)+model_abss25_response(3,:))/2,'linewidth',2,'color',model_color*1.2);
legend([h3;h4;h1;h2],'Meas. 1.7 mab','Meas. 0.15 mab','Model 1.6 mab','Model 0.15 mab')
%ylim([0 .05])
ylabel('Concentration (kg/m^3)')
title('2.5 MHz ABSS')
legend([h3;h1],'Data','Model')

pfn = sprintf('abss_data_model_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% ABSS skill
% index in elev
iizhi = find(elev>=1.9,1,'first');
fprintf(1,'in floc_plot, iiz = %d and elev(iiz) = %6.2f\n',iizhi,elev(iizhi))
iizlo = find(elev>=0.2,1,'first');
fprintf(1,'in floc_plot, iiz = %d and elev(iiz) = %6.2f\n',iizlo,elev(iizlo))

fprintf(1,'Skill for modeled and measured 2.5 MHz ABSS at z=%6.2f\n',elev(iizhi))
imodel = interp1(s2d*ocean_time,model_abss25_response(iizhi,:),ydd,'nearest');
s1 = skill(imodel,NX(:,24))
[rmsd_star,bias,r]=target_diagram(imodel,NX(:,24));
fprintf(fid,'%6.2f, %6.2f, %6.3f, %6.2f, %6.2f, %6.2f, ',s1.willmott,s1.brier,s1.RMS,rmsd_star,bias,r)

figure(9); clf
plot(NX(:,24),imodel,'ok')
hold on
fprintf(1,'Skill for modeled and measured 2.5 MHz ABSS at z=%6.2f\n',elev(iizlo))
imodel = interp1(s2d*ocean_time,model_abss25_response(iizlo,:),ydd,'nearest');
s2 = skill(imodel,NX(:,27))
[rmsd_star,bias,r]=target_diagram(imodel,NX(:,27));
fprintf(fid,'%6.2f, %6.2f, %6.3f, %6.2f, %6.2f, %6.2f, ',s2.willmott,s2.brier,s2.RMS,rmsd_star,bias,r)
plot(NX(:,24),imodel,'or')
axis([0 .01 0 .01])
xlabel('Data')
ylabel('Model')
axis('square')
pfn = sprintf('abss_data_model_scatter_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end

%% ac9 trans

% NX(3) is acn9 at 2 mab, NX(4) is ac9 at .9 mab
% convert from 1/m to kg/m3  as done in ws_fit: allc = (allc*2.7)-1.3;
figure(10); clf
h3=plot(ydd,.01*(2.7*NX(:,3)-1.3),'linewidth',2);
hold on
set(h3,'color',[.2 .2 .2])
h4=plot(ydd,.01*(2.7*NX(:,4)-1.3),'linewidth',2);
set(h4,'color',[.4 .4 .4])
h1=plot( (18./24.)+s2d*ocean_time,model_ac9_response(10,:),'linewidth',2);
h2=plot( (18./24.)+s2d*ocean_time,model_ac9_response(5,:),'linewidth',2);
title('ac9 SSC Estimates')
ylabel('Concentration (kg/m^3)')
ylim([ 0 .30 ])
legend([h3;h4;h1;h2],'Data, 2 mab','Data, 1 mab','Model 0.9 mab','Model 2.1 mab') 
pfn = sprintf('ac9_data_model_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% Compare ac9 ws
figure(11); clf
h6=plot( (18./24.)+s2d*ocean_time,-1e3*[pfac9.p]'.*(0.41*ustrc),'linewidth',2);
hold on 
ok = find(UMEattn650_r22h>0.5);
h2=plot(ydd(ok),1000*UMEattn650_ws2h(ok),'linewidth',2);
title('ac9 Settling Velocity')
pfn = sprintf('ac9_ws_data_model_ts_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end

%% align ac9 model and data, compute skill
% index in elev
iizhi = find(elev>=2.0,1,'first');
fprintf(1,'in floc_plot, iiz = %d and elev(iiz) = %6.2f\n',iizhi,elev(iizhi))
iizlo = find(elev>=0.9,1,'first');
fprintf(1,'in floc_plot, iiz = %d and elev(iiz) = %6.2f\n',iizlo,elev(iizlo))

fprintf(1,'Skill for modeled and measured ac9 at z=%6.2f\n',elev(iizhi))
imodel = interp1(s2d*ocean_time,model_ac9_response(iizhi,:),ydd,'nearest');
s1 = skill(imodel,2.7*NX(:,3)-1.3)
[rmsd_star,bias,r]=target_diagram(imodel,2.7*NX(:,3)-1.3);
fprintf(fid,'%6.2f, %6.2f, %6.3f, %6.2f, %6.2f, %6.2f, ',s1.willmott,s1.brier,s1.RMS,rmsd_star,bias,r)

figure(12); clf
h1=plot(2.7*NX(:,3)-1.3,imodel,'ok')
hold on

fprintf(1,'Skill for modeled and measured ac9 at z=%6.2f\n',elev(iizlo))
imodel = interp1(s2d*ocean_time,model_ac9_response(iizlo,:),ydd,'nearest');
s2 = skill(imodel,2.7*NX(:,4)-1.3)
[rmsd_star,bias,r]=target_diagram(imodel,2.7*NX(:,4)-1.3);
fprintf(fid,'%6.2f, %6.2f, %6.3f, %6.2f, %6.2f, %6.2f, ',s2.willmott,s2.brier,s2.RMS,rmsd_star,bias,r)
h2=plot(2.7*NX(:,4)-1.3,imodel,'or')
%axis([0 .01 0 .01])
xlabel('Data')
ylabel('Model')
axis('square')
pfn = sprintf('ac9_data_model_scatter_%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end



