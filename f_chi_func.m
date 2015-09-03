function [ffi,chisvo,x] = f_chi_func( a, f, rhof )
% function [ffi,chisvo,x] = f_chi_func( a, f, rhof )
% Compute intrinsic form and scattering functions
%
% Input:
%   a - particle radius (m)
%   f - acoustic frequency (Hz)
%   rhof - particle density (kg/m3) [optional - defaults to rhos]
%
% Returns:
%   ffi - intrinsic form function
%   x   - nondimensional particle size

% Default to solid particles
global mvco
% These physcial parameters could be passed as arguments or global variables
% For comparison with Thorne:
co=1480;   % speed of sound in water
rhow=1000; % density of water
ck=5550;   % speed of sound in particle
rhos=2600; % density of particle
v=1e-6;    % kinematic  
% For MVCO:
if(mvco)
   co=1500.;   % speed of sound in water
   rho0=1025.; % density of water
   ck=5800;    % speed of sound in quartz
   rhos=2650;  % density of quartz  
end

x=2*pi*a*f/co;
Co=1/(rhow*co^2);
Ck=1/(rhos*ck^2); %compressibility
gk=rhos/rhow;
hk=ck/co;

% effective density
if(exist('rhof','var')~=1), rhof = rhos;, end
rhoex = rhof-rhow;

% gamma_o  zeta_o values
beta1=1.02; % gamma_o
beta2=1.02; % zeta_o

% maximum physical value for effective density
rhoex = min(rhoex, rhos-rhow);
% not letting effective density go below zero
rhoex = max( 0., rhoex );
gam=1.+rhoex/rhow;
gam = max( gam, beta1 )
%zz=find(gam<beta1); gam(zz)=beta1;
phi=(gk-gam)/(gk-1);

% Wood's expression - Eqn 9
vw=((phi*rhow+(1-phi)*rhos).*(phi*Co+(1-phi)*Ck)).^(-0.5);
h=vw/co;
zz=find(h<beta2); h(zz)=beta2;

% form function f_fi
% unnumbered Clay & Medwin eqn on p. 65
e=(gam.*h.^2);
kf=2*abs((e-1)./(3*e)+(gam-1)./(2*gam+1));
alpha1=1.2; % called epsilon_1 on p. 85
ffi=(kf.*x.^2)./(1+alpha1*x.^2); % Eqn. 6a

% second Clay & Medwin eqn. on p. 85
kalfa=2*(((e-1)./(3*e)).^2+(1/3)*((gam-1)./(2*gam+1)).^2);

% chi_h scattering; Eqn. 6b
chishc=(kalfa.*x.^4)./(1-1.0*x+1.5*x.^2+kalfa.*x.^4);
%subplot(2,1,2), loglog(x,chishc,':k'), hold on

%viscous atten
% chi_hs scattering

% rho=g; %rhok/rho0; 


% da=c*dx/(2*pi*f); 
% a=a(1:nz1-1); 
w=2*pi*f;
k=2*pi*f/co;
beta=sqrt(w/(2*v));
theta=0.5*(1+(9./(2*beta*a))); % after eqn 3c
tau=(9./(4*beta*a)).*(1+(1./(beta*a)));
e1=((k*(gam-1).^2)/2);
e2=tau./(tau.^2+(gam+theta).^2);

%converted to chi_hv values
e12=e1.*e2;
chiv=(4*a/3).*e12;

%combine scat and visc atten chi_h
chisvo=chishc+chiv;