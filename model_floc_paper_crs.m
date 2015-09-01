%n:\mat5work\sizescatter\floc\model_floc_paper.m 26/06/14

%This software is in support of the paper
% 'Modelling acoustic scattering by suspended flocculating sediments'
%by Peter D. Thorne1, Iain T. MacDonald and Christopher E. Vincent in Continental Shelf Research

% The program is for educational use and not for comercial purposes.
% The code is unsupported.


%1 The form function is calculated for the primary, intermediate and fully formed flocs as a function of particle size
%2 The average form function is calculated for the primary, intermediate and fully formed flocs as a function of particle size
%3 The normalise total scattering and viscous cross-section is calculated for the primary, intermediate and fully formed flocs as a function of particle size
%4 The average normalise total scattering and viscous cross-section is calculated for the primary, intermediate and fully formed flocs as a function of particle size
%5 Plots of the results are given
% close all
% clear all
f = 3e6;

nz1=1e5;
a=logspace(-10,-2,nz1);

rhow=1000;co=1480;rhos=2650;
x=2*pi*a*f/co;
dx=diff(x); x=x(1:nz1-1); a=a(1:nz1-1);
nz=100;
aa=logspace(-9,-3,nz);
xx=2*pi*aa*f/co;


% function [f,chi,x]=intrinsic_fchi( f, a, rhof, rhos, rhow )
% defining primary parameters

ck=5550; 
Co=1/(rhow*co^2); % compressibility
Ck=1/(rhos*ck^2);
gk=rhos/rhow;
hk=ck/co;


% Follwing formulations are for form of rhoex( a );
% rhoex=1e-3./(a.^1); %choosen for illustration
% limit to realistic values
% rhoex(rhoex>(rhos-rhow))=rhos-rhow;

% effective density
rhoex = rhos-rhow;
rhoex(rhoex<0)=0;

% gamma_o  zeta_o values
beta1=1.02; % gamma_o
beta2=1.02; %zeta_o

g=1+rhoex/(rhow);
g(g<beta1)=beta1;
phi=(gk-g)/(gk-1);

% Calculate h (zeta in Eqn. 9) using Wood expression
vw=((phi*rhow+(1-phi)*rhos).*(phi*Co+(1-phi)*Ck)).^(-0.5);
h=vw/co;
%zz=find(h<beta2); h(zz)=beta2;
h(h<beta2)=beta2;
e=(g.*h.^2);

% **************************   form function   ************************
% form function f_h

kf=2*abs((e-1)./(3*e)+(g-1)./(2*g+1)); % Eqn 
alpha1=1.2;
fmj=(kf.*x.^2)./(1+alpha1*x.^2); % Eqn. 6a

% avg form function for lognormal distribution f_ho
dsigma=0.5;
%Lognormal
for jj=1:nz;
   xo=xx(jj);
   sigma=dsigma*xo;
   mu=log((xo.^2)./sqrt(xo.^2+sigma.^2));
   sigman=sqrt(log((sigma/xo).^2+1));
   plognorm=(1./(x*sigman*sqrt(2*pi))).*exp(-((log(x)-mu).^2)/(2*sigman.^2));
   axlng(jj)=sum(x.*plognorm.*dx);
   ax3=sum((x.^3).*plognorm.*dx);
   axf=sum(((x.^2).*(fmj.^2).*dx).*plognorm);
   aflng(jj)=sqrt((axlng(jj).*axf)./ax3);
end

mn=15; % font size
figure;,orient tall
subplot(2,1,1),
loglog(x,fmj,'k'), hold on
xlabel('x_o','fontsize',mn), ylabel('f_{ho}','fontsize',mn)
axis([1e-4 1e1 1e-7 1e-0])
set(gca,'fontsize',mn)
text(5,0.4,'a','fontsize',mn)
loglog(axlng,aflng,'--k')
hh=legend(' \delta=0.0','\delta=0.5');
set(hh,'fontsize',12)

%% **************************   atten   ************************
% chi_h scattering
kalfa=2*(((e-1)./(3*e)).^2+(1/3)*((g-1)./(2*g+1)).^2);
chishc=(kalfa.*x.^4)./(1-1.0*x+1.5*x.^2+kalfa.*x.^4);
%subplot(2,1,2), loglog(x,chishc,':k'), hold on

%viscous atten
% chi_hs scattering
rho=g; %rhok/rho0
nu=1e-6; % kinematic viscosity

da=co*dx/(2*pi*f);
a=a(1:nz1-1);
w=2*pi*f;       % angular frequency
k=2*pi*f/co;    % wave number
beta=sqrt(w/(2*nu));
delta=0.5*(1+(9./(2*beta*a)));
s=(9./(4*beta*a)).*(1+(1./(beta*a)));
e1=((k*(rho-1).^2)/2);
e2=s./(s.^2+(rho+delta).^2);

%converted to chi_hv values
e12=e1.*e2;
chiv=(4*a/3).*e12;

%combine scat and visc atten chi_h
chisvo=chishc+chiv;

% ave scattering for lognormal distribution chi_ho
% use dsigma from above
%Lognormal
for jj=1:nz;
   xo=xx(jj);
   sigma=dsigma*xo;
   mu=log((xo.^2)./sqrt(xo.^2+sigma.^2));
   sigman=sqrt(log((sigma/xo).^2+1));
   plognorm=(1./(x*sigman*sqrt(2*pi))).*exp(-((log(x)-mu).^2)/(2*sigman.^2));
   axlng(jj)=sum(x.*plognorm.*dx);
   ax3=sum((x.^3).*plognorm.*dx);
   axchi=sum((x.^2).*chisvo.*plognorm.*dx);
   achilngo(jj)=(axlng(jj).*axchi)./ax3;
end

subplot(2,1,2)
loglog(x,chisvo,'k'), hold on
axis([1e-4 1e1 1e-7 1e-0])
xlabel('x_o','fontsize',mn), ylabel('\chi_{ho}','fontsize',mn)
set(gca,'fontsize',mn)
%title('Normalised total scattering cross-section  - attenuation','fontsize',mn)
text(5,0.4,'b','fontsize',mn)
loglog(axlng,achilngo,'--k')

