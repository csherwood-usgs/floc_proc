clear

%% case number
cas='38'


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
%rho=1025;
rho=nanmean(nanmean(nc{'rho'}(:,:,3,3)))+1000;
bstru=nc{'bustrcwmax'}(:,3,3);
bstrv=nc{'bvstrcwmax'}(:,3,3);
bstr=abs(complex(bstru,bstrv));
bstruc=nc{'bustrc'}(:,3,3);
bstrvc=nc{'bvstrc'}(:,3,3);
bstrc=abs(complex(bstruc,bstrvc));
ubar=nc{'ubar'}(:,3,3);
vbar=nc{'vbar'}(:,3,3);
speed=abs(complex(ubar,vbar));
Zo_app=nc{'Zo_app'}(:,3,3);
%%
figure(3);clf
line(tim,bstr)
line(ustar_av.dn,ustar_av.ustrr.*ustar_av.ustrr*rho,'Color','r')
line(tim,bstrc,'Color','k','LineStyle','--')
line(ustar_av.dn,ustar_av.ustrc.*ustar_av.ustrc*rho,'Color','k')
axtt;legend('w-c model ','w-c obs.','curr. mod','curr. obs',0)
datetick('keeplimits','keepticks')
title('stress')
%%
figure(2);clf
line(tim,bstrc,'Color','k','LineStyle','--')
line(ustar_av.dn,ustar_av.ustrc.*ustar_av.ustrc*rho,'Color','r')
axtt;legend('current model ','current obs.',0)
%datetick('keeplimits','keepticks')
title('stress')
%%
figure(4);clf
line(tim,Zo_app)
line(ustar_av.dn,ustar_av.zoa,'Color','r')
axtt;legend('zoa model ','zoa obs.',0)
datetick('keeplimits','keepticks')
title('apparent roughness')

%%
a=(ustar_av.ustrc.*ustar_av.ustrc*rho);
a(1)=a(2);
b=bstrc;
cor=lagcor(a,b,6)    
cor2=lagcor(b,a,6)
ai=find(cor2==max(cor2))
%%

figure(22);clf
line(tim-tim(3)+tim(1),bstrc,'Color','k','LineStyle','--')
line(ustar_av.dn,ustar_av.ustrc.*ustar_av.ustrc*rho,'Color','r')
axtt;legend('current model ','current obs.',0)
%datetick('keeplimits','keepticks')
title('current stress')



