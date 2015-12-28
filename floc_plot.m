% floc_plot - makes plots after floc_proc
set_colors
s2d = 1. /(3600.*24.);
% index in elev
iiz = find(elev>=0.5,1,'first');
fprintf(1,'in floc_plot, iiz = %d and elev(iiz) = %f\n',iiz,elev(iiz))
fprintf(fid,'%d, %d, %d, %3.1f, ',nz,NCS,iiz,elev(iiz))
% text string for plot labels
ets = sprintf(' at z = %4.1f mab',elev(iiz));
%% floc conc plots - pcolor version = muds, wst, fdiamt
figure(9); clf
subplot(311)
pcolorjw( s2d*tz, h+z_w, log10(muds+eps))
caxis([-2.5 1.5])
set(gca,'xticklabel','','fontsize',14)
colorbar
title('log_{10} Total Floc Concentration (kg/m^3)')

% fraction-weighted ws
wst = squeeze(sum(repmat(ws(1:NCS),1,nz,nt).*m)./sum(m));
subplot(312)
pcolorjw( s2d*tz, h+z_w, 1e3*wst)
set(gca,'xticklabel','','fontsize',14)
colorbar
title('Floc Settling Velocity (mm/s)')
ylabel('Elevation (m)')

% fraction-weighted size
subplot(313)
fdiamt = squeeze(sum(repmat(fdiam(1:NCS),1,nz,nt).*m)./sum(m));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamt)
colorbar
title('Floc Diameter (\mum)')
set(gca,'fontsize',14)
xlabel('Days','fontsize',16)
pfn=sprintf('floc_conc_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end

%% sand conc plots - snds, wssdt, fdiamsdt
figure(10); clf
subplot(311)
pcolorjw( s2d*tz, h+z_w, log10(snds+eps))
set(gca,'xticklabel','','fontsize',14)
caxis([-2.5 1.5])
colorbar
title('log_{10} Total Sand Concentration (kg/m^3)')

% fraction-weighted ws
wssdt = squeeze(sum(repmat(ws(NNN+NCS+1:NCS+NNN+NND),1,nz,nt).*snd)./sum(snd));
subplot(312)
pcolorjw( s2d*tz, h+z_w, 1e3*wssdt)
set(gca,'xticklabel','','fontsize',14)
colorbar
ylabel('Elevation (m)','fontsize',16)
title('Sand Settling Velocity (mm/s)')

% fraction-weighted size
subplot(313)
fdiamsdt = squeeze(sum(repmat(fdiam(NNN+NCS+1:NCS+NNN+NND),1,nz,nt).*snd)./sum(snd));
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamsdt)
colorbar
set(gca,'fontsize',14)
xlabel('Days','fontsize',16)
title('Sand Diameter (\mum)')
pfn=sprintf('sand_conc_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end

%% fraction-weighted size of both sand and flocs - fdiamall
figure(11); clf
subplot(211)
h1=plot(s2d*ocean_time, 1e6*fdiamsdt(iiz,:) );
hold on
h2=plot(s2d*ocean_time, 1e6*fdiamt(iiz,:) );
h3=plot(s2d*ocean_time, 1e6*fdiamall(iiz,:) );
legend([h1;h2;h3],'Sand','Flocs','Combined')
ylabel('Mass-weighted Mean Diameter (\mum)')
subplot(212)
pcolorjw( s2d*tz, h+z_w, 1e6*fdiamt)
colorbar
title(['Floc + Sand Diameter (\mum)',ets])
set(gca,'fontsize',14)
xlabel('Days','fontsize',16)
pfn=sprintf('floc+sand_size_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%fprintf(fid,'Mean conc and size at %f mab:\n  flocs: %f %5.1f\n   sand: %f %5.1f\n',...
fprintf(fid,'%7.4f, %5.1f, %7.4f, %5.1f\n',...
   nanmean(muds(iiz,:)),nanmean(1e6*fdiamt(iiz,:)),...
   nanmean(snds(iiz,:)),nanmean(1e6*fdiamsdt(iiz,:)))

%% line plots of concs and diam
figure(12); clf
h1=plot(s2d*ocean_time, snds(iiz,:),'linewidth',2,'color',sand_color);
hold on
h2=plot(s2d*ocean_time, muds(iiz,:),'linewidth',2,'color',floc_color);
h3=plot(s2d*ocean_time, (snds(iiz,:)+muds(iiz,:)),'linewidth',2,'color',comb_color);
% iiz=9
% h4=plot(s2d*ocean_time, snds(iiz,:),'--','linewidth',2,'color',sand_color)
% hold on
% h5=plot(s2d*ocean_time, muds(iiz,:),'--','linewidth',2,'color',floc_color)
% h6=plot(s2d*ocean_time, (snds(iiz,:)+muds(iiz,:)),'--','linewidth',2,'color',comb_color);
title(['Modeled Concentration of Flocs and Sand',ets])
legend([h1;h2;h3],'Sand','Flocs','Combined')

%% line plots of combined conc and diam
figure(13); clf
h1=plot(s2d*ocean_time, vcomb1(iiz,:),'linewidth',2,'color',abss1c_color);
hold on
h2=plot(s2d*ocean_time, vcomb4(iiz,:),'linewidth',2,'color',abss3c_color);
% h2=plot(s2d*ocean_time, vcomb25(iiz,:),'linewidth',2,'color',abss2c_color)
h3=plot(s2d*ocean_time, (r_ac9f(3,:)+r_ac9s(iiz,:)+eps),'linewidth',2,'color',ac9c_color);
h4=plot(s2d*ocean_time, (r_lisstf(3,:)+r_lissts(iiz,:)+eps),'linewidth',2,'color',lisstc_color);
title(['Combined Response to Flocs and Sand', ets])
legend([h1;h2;h3;h4],'1 MHz ABSS','4 MHz ABSS','ac9','LISST')
%% Plot acoustic response
figure(14); clf
subplot(311)
pcolorjw( s2d*tz, h+z_w, v25)
colorbar
title('Flocs Acoustic Response at 2.5 MHz')
set(gca,'xticklabel','','fontsize',14)
%xlabel('Days')
%ylabel('Elevation (m)')
xlim([0 32.8])
caxis([ 0 .5])
ylim([0 3])
subplot(312)
pcolorjw( s2d*tz, h+z_w, vsn25)
colorbar
title('Sand Acoustic Response')
set(gca,'xticklabel','','fontsize',14)
%xlabel('Days')
ylabel('Elevation (m)','fontsize',16)
xlim([0 32.8])
caxis([ 0 .5])
ylim([0 3])

subplot(313)
pcolorjw( s2d*tz, h+z_w, vcomb25)
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

%% Plot optical responses - pcolor version
figure(15); clf
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

pfn=sprintf('pcol_both_optical_response_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% line version
figure(16); clf;
h1=plot(s2d*ocean_time, r_ac9s(iiz,:),'linewidth',2,'color',ac9s_color);
hold on
h2=plot(s2d*ocean_time, r_ac9f(iiz,:),'linewidth',2,'color',ac9f_color);
h3=plot(s2d*ocean_time, r_lissts(iiz,:),'linewidth',2,'color',lissts_color);
h4=plot(s2d*ocean_time, r_lisstf(iiz,:),'linewidth',2,'color',lisstf_color);
title(['Optical Response to Flocs and Sand',ets])
ylabel('Concentration (kg/m^3)')
legend([h1;h2;h3;h4],'ac9 sand','ac9 flocs','LISST sand','LISST flocs');
pfn=sprintf('line_both_optical_response_run%02d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end
%% time series of sand, flocs and combined masses
figure(17); clf
h1=plot(s2d*ocean_time, vcomb1(iiz,:),'linewidth',2,'color',abss1c_color);
hold on
%h2=plot(s2d*ocean_time, vcomb4(iiz,:),'linewidth',2,'color',abss3c_color);
h3=plot(s2d*ocean_time, (r_ac9f(3,:)+r_ac9s(iiz,:)+eps),'linewidth',2,'color',ac9c_color);
h4=plot(s2d*ocean_time, (r_lisstf(3,:)+r_lissts(iiz,:)+eps),'linewidth',2,'color',lisstc_color);
title(['Combined Response to Flocs and Sand',ets])
legend([h1;h3;h4],'1 MHz ABSS','ac9','LISST')
%% time series of both combined responses
% TODO Edit this to plot estimated size
% figure(10); clf
% h1=plot(s2d*ocean_time, vcomb1(iiz,:),'linewidth',2,'color',abss1c_color);
% hold on
% h2=plot(s2d*ocean_time, vcomb4(iiz,:),'linewidth',2,'color',abss3c_color);
% h3=plot(s2d*ocean_time, (r_ac9f(3,:)+r_ac9s(iiz,:)+eps),'linewidth',2,'color',ac9c_color);
% h4=plot(s2d*ocean_time, (r_lisstf(3,:)+r_lissts(iiz,:)+eps),'linewidth',2,'color',lisstc_color);
% title(['Combined Response to Flocs and Sand',ets])
% legend([h1;h2;h3;h4],'1 MHz ABSS','4 MHz ABSS','ac9','LISST')
%%
load cmap_plusminus
figure(18)
pcolorjw( s2d*tz, h+z_w, vcomb25-vsn25)
colorbar
colormap(cmap_plusminus)
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
figure(19); clf
h1=plot(tsed(:),vcomb25(:),'.','color',abss2c_color);
set(h1,'markersize',14);
hold on
h2=plot(muds(:),v25(:),'.','color',abss2f_color);
set(h2,'markersize',14);
h3=plot(snds(:),vsn25(:),'.','color',abss2s_color);
set(h3,'markersize',14);
xlabel('Mass Concentration (kg/m^3)','fontsize',16)
ylabel('Acoustic Response','fontsize',16)
h4=legend([h3;h2;h1],'Sand','Flocs','Combined');
set(h4,'fontsize',14)
pfn = sprintf('acoustic_response_scatter%2d.png',cas)
if(iplot),print('-dpng','-r300',pfn); end