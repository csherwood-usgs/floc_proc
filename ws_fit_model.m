% ws_fit_model - fit Rouse profiles to
iplot=1;
za=0.1
cc=r_ac9f+r_ac9s+eps;
for ii = 1:nt
   c = cc(1:11,ii);
   z=h+z_r(1:11,ii);
   pfoptics(ii) = pfit( c, z,1, za);
end


for ii = 1:nt
   c = v(1:11,ii);
   z=h+z_r(1:11,ii);
   pfacoustics(ii) = pfit( c, z,1, za);
end

for ii = 1:nt
   c = vcomb(1:11,ii);
   z=h+z_r(1:11,ii);
   pfacoustics_comb(ii) = pfit( c, z,1, za);
end

%% compare with observations
rho0 = ncread(url,'rho0')
tauc = (squeeze(ncread(url,'bustrc',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrc',[i j 1],[1 1 Inf])).^2).^(.5);
taucwm = (squeeze(ncread(url,'bustrcwmax',[i j 1],[1 1 Inf])).^2+squeeze(ncread(url,'bvstrcwmax',[i j 1],[1 1 Inf])).^2).^(.5);
%%
ustar = sqrt( tauc/rho0 );
wsoptics = 0.41*ustar.*[pfoptics.p]';
wsacoustics = 0.41*ustar.*[pfacoustics.p]';
wsacoustics_comb = 0.41*ustar.*[pfacoustics_comb.p]';
%%
tt = s2d*tz(1,:)';
figure(1);clf
h1=plot(tt,-1e3*wsacoustics_comb,'linewidth',2);
set(h1,'color',[.6 .6 .6]);
set(gca,'fontsize',14)
hold on
h3=plot(tt,-1e3*wsoptics,'linewidth',2);
set(h3,'color',[.8 .4 .4]);
h2=plot(tt,-1e3*wsacoustics,'linewidth',2);
set(h2,'color',[.6 .6 1]);

h4=legend([h1;h2;h3],'Combined Acoustics','Optics','Acoustics (sand)')
xlabel('Days','fontsize',16)
ylabel('Apparent {\itw_s} (mm/s)','fontsize',16);
xlim([0 32.8])
pfn=sprintf('model_ws_timeseries_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% look for size inversions
dfd2 = 1e6*(fdiamt(9,:)-fdiamt(2,:));
dfdall = 1e6*(fdiamall(9,:)-fdiamall(2,:));
figure(2); clf
plot([0 tt(end)],[0 0],'--k')
hold on
h1=plot(tt,dfd2,'linewidth',2,'color',[.6 .6 .6]);
hold on
h2=plot(tt,dfdall,'linewidth',2,'color',[.3 .3 .3]);
set(gca,'fontsize',14)
ylabel('{\itD}_{1.9} - {\itD}_{0.6} ({\mu}m)','fontsize',16)
xlabel('Days','fontsize',16)
xlim([0 32.8])
pfn=sprintf('model_delta_D_timeseries_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end

figure(3); clf
plot([0 s2d*tz(end)],[0 0],'--k')
hold on
h1=plot(taucwm,dfd2,'.')
set(gca,'fontsize',14)
set(h1,'color',[.2 .2 .2],'markersize',12)
ylabel('{\itD}_{1.9} - {\itD}_{0.6} ({\mu}m)','fontsize',16)
xlabel('Wave-current Combined Stress (Pa)')
xlim([0 5])
pfn=sprintf('model_delta_D_stress_%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end