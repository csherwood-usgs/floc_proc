clear

%% case number
cas='52'


%%

load('ustar_av')

url=['http://geoport.whoi.edu/thredds/dodsC/clay/usgs/users/aretxabaleta/MVCO/ocean_his_',cas,'.nc'];

nc=ncgeodataset(url);

tim=nc{'ocean_time'}(:)/3600/24+ustar_av.dn(1);
nt=length(tim);

zeta = squeeze( nc{'zeta'}(:,3,3) );
h        = nc{'h'}(3,3);
Cs_r = nc{'Cs_r'}(:) ;    Cs_w = nc{'Cs_w'}(:) ;
nz=length(Cs_r);
s_w = nc{'s_w'}(:) ;    s_rho = nc{'s_rho'}(:) ;    hc = nc{'hc'}(:) ;
for  k = 1:nz+1
    zw(:,k) = zeta.*(1+s_w(k)*hc./h-hc*Cs_w(k)./h+Cs_w(k))+(h - hc)*Cs_w(k);
end
for  k = 1:nz
  z_r(:,k) = zeta.*(1+s_rho(k)*hc./h-hc*Cs_r(k)./h+Cs_r(k))+hc*s_rho(k)+(h - hc)*Cs_r(k);
end
dz = diff(zw,1,2);

rho=nanmean(nanmean(nc{'rho'}(:,:,3,3)))+1000;
bstru=nc{'bustrcwmax'}(:,3,3);
bstrv=nc{'bvstrcwmax'}(:,3,3);
bstr=abs(complex(bstru,bstrv));
%%
vonKar=0.41;
gls=nc{'gls'}(:,:,3,3);
tke=nc{'tke'}(:,:,3,3);
gls_av(:,1)=gls(:,1);
tke_av(:,1)=tke(:,1);
for k=2:nz
    gls_av(:,k)=.5*gls(:,k)+.5*gls(:,k+1);
    tke_av(:,k)=.5*tke(:,k)+.5*tke(:,k+1);
end
% for k=1:nz
%     gls_av(:,k)=.5*gls(:,k)+.5*gls(:,k+1);
%     tke_av(:,k)=.5*tke(:,k)+.5*tke(:,k+1);
% end
gls_p = -1  ;%                        ! k-omega
gls_m = 0.5;
gls_n = -1.0;
gls_cmu0 = 0.5477;
exp1 = 3.0+gls_p/gls_n;
exp2 = 1.5+gls_m/gls_n;
exp3 = -1.0/gls_n;
diss0 = (gls_cmu0^exp1).*(tke_av.^exp2).*(gls_av.^exp3);
diss=diss0;
diss(:,1)=(bstr.^1.5)/(vonKar*0.9);

nu0=1.5e-6;
G0=sqrt(diss0/nu0);
G=sqrt(diss/nu0);
%%
figure(1);clf
pcolorjw(tim*ones(1,nz),z_r,G)
caxis([0,20])
datetick('keepticks','keeplimits')
colorbar

%%
figure(2);clf
set(gcf,'PaperPosition',[.5,.5,10,8]);wysiwyg
subplot 221
pcolorjw(tim*ones(1,nz),z_r,G0)
caxis([0,4])
datetick('x','mm/dd','keepticks','keeplimits')
colorbar
title('G from GLS')
subplot 222
pcolorjw(tim*ones(1,nz),z_r,G)
caxis([0,4])
datetick('x','mm/dd','keepticks','keeplimits')
colorbar
title('enhanced G ')
subplot 223
pcolorjw(tim*ones(1,nz),z_r,G0)
set(gca,'YLim',[-11.9,-10])
caxis([0,100])
datetick('x','mm/dd','keepticks','keeplimits')
colorbar
title('G from GLS zoom ')
subplot 224
pcolorjw(tim*ones(1,nz),z_r,G)
set(gca,'YLim',[-11.9,-10])
caxis([0,2000])
datetick('x','mm/dd','keepticks','keeplimits')
colorbar
title('enhanced G zoom ')
eval(['print -dpng -painters G_',cas,'.png'])
