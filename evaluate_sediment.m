% evaluate_sediment.m - Script to check sediment conservation
% requires nctoolbox
clear
%% case number
cas='78'
%%
load('ustar_av')

url=['http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/aretxabaleta/MVCO/ocean_his_',cas,'.nc'];

nc=ncgeodataset(url);

tim=nc{'ocean_time'}(:)/3600/24+ustar_av.dn(1);
nt=length(tim);

zeta = squeeze( nc{'zeta'}(:) );
h        = nc{'h'}(:);
Cs_r = nc{'Cs_r'}(:) ;    Cs_w = nc{'Cs_w'}(:) ;
nz=length(Cs_r);
s_w = nc{'s_w'}(:) ;    s_rho = nc{'s_rho'}(:) ;    hc = nc{'hc'}(:) ;
[nx,ny]=size(h);
%%
for i=1:nx
   for j=1:ny
      for  k = 1:nz+1
         zw(:,k,i,j) = squeeze(zeta(:,i,j)).*(1+s_w(k)*hc./h(i,j)-hc*Cs_w(k)./h(i,j)...
            +Cs_w(k))+(h(i,j) - hc)*Cs_w(k);
      end
      for  k = 1:nz
         z_r(:,k,i,j) = squeeze(zeta(:,i,j)).*(1+s_rho(k)*hc./h(i,j)-hc*Cs_r(k)./h(i,j)...
            +Cs_r(k))+hc*s_rho(k)+(h(i,j) - hc)*Cs_r(k);
      end
   end
end
dz = diff(zw,1,2);
%% get sediment

nsed=15;

% Get cohesive fraction and mass from netcdf files
for iSed = 1:nsed
   iStr = sprintf('%02d',iSed);
   mud(iSed,:,:,:,:)     = nc{['mud_' iStr]}(:);
   mudmass(iSed,:,:,:,:) = nc{['mudmass_' iStr]}(:);
   mudfrac(iSed,:,:,:,:) = nc{['mudfrac_' iStr]}(:);
end


%% Calculate total suspended and bed mass
nbed=size(mudfrac,3);
% Initialize variables for combined arrays
susmass   = zeros(nsed,nt,nx,ny);       % Depth-integrated suspended mass
susmass_z = zeros(nsed,nt,nz,nx,ny);    % Depth-variable suspended mass
frac      = zeros(nsed,nt,nbed,nx,ny);  % Bed sediment fraction
sedmass   = zeros(nsed,nt,nbed,nx,ny);  % Bed sediment mass
sedmass_t = zeros(nsed,nt,nx,ny);       % Depth-integrated bed sediment mass

% Calculate totals in bed per class
mudmass_t  = sum(mudmass,3);    % total mass of each mud class in bed
%%
for iSed=1:nsed
   sdz = squeeze(mud(iSed,:,:,:,:)).*dz;
   
   susmass(iSed,:,:,:)     = squeeze(sum(sdz,2));
   susmass_z(iSed,:,:,:,:) = mud(iSed,:,:,:,:);
   
   sedmass(iSed,:,:,:,:)   = mudmass(iSed,:,:,:,:);
   frac(iSed,:,:,:,:)      = mudfrac(iSed,:,:,:,:);
   sedmass_t(iSed,:,:,:)   = squeeze(mudmass_t(iSed,:,:,:,:));
end

susmasstot   = squeeze(sum(susmass,1));      % total mass in suspension
susmasstot_z = squeeze(sum(susmass_z,1));    % total mass in suspension per level
sedmasstot   = squeeze(sum(sedmass_t,1));    % total mass in bed

%%
total_mass = sedmass_t+susmass;          % total sediment mass (should be conserved)

totalsed=squeeze(sum(total_mass,1));
dum=squeeze(sum(totalsed,3));
toto=squeeze(sum(dum,2));
%%

figure(3);clf
line(tim,toto)
datetick('keeplimits','keepticks')
eval(['print -dpng total_sed_',cas,'.png'])





