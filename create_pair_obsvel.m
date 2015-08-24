clear

load('ustar_av')

%% use stress right dx with just currents
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);

dx = (x_rho(1,2)-x_rho(1,1))/10;
dpdx = rho*ustar_av.ustrc.*ustar_av.ustrc;
dpdx(1)=dpdx(2);
dpdx=dpdx.*sign(ustar_av.ur);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_sign_rot_velobs_curr_stress.nc', 'Pair', 'Pair')





%%
figure;
line(ocean_time,dpdx)

%% use stress right dx with just currents
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);

dx = (x_rho(1,2)-x_rho(1,1))/5;
dpdx = rho*ustar_av.ustrc.*ustar_av.ustrc;
dpdx(1)=dpdx(2);
dpdx=dpdx.*sign(ustar_av.ur);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_sign_rot_2_curr_stress.nc', 'Pair', 'Pair')

%% use stress right dx with just currents shifted 2 hours
ocean_time=ustar_av.dn-ustar_av.dn(3);
rho=1025;
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);

dx = (x_rho(1,2)-x_rho(1,1))/5;
dpdx = rho*ustar_av.ustrc.*ustar_av.ustrc;
dpdx(1)=dpdx(2);
dpdx=dpdx.*sign(ustar_av.ur);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_sign_rot_2_shift_curr_stress.nc', 'Pair', 'Pair')



%% use stress right dx with just currents
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);

dx = (x_rho(1,2)-x_rho(1,1));
dpdx = rho*ustar_av.ustrc.*ustar_av.ustrc;
dpdx(1)=dpdx(2);
dpdx=dpdx.*sign(ustar_av.ur);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_sign_rot_10_curr_stress.nc', 'Pair', 'Pair')



return
%% use stress right dx with just currents
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);

dx = (x_rho(1,2)-x_rho(1,1))/10;
dpdx = rho*ustar_av.ustrc.*ustar_av.ustrc;
dpdx(1)=dpdx(2);
dpdx=dpdx.*sign(ustar_av.u);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_sign_velobs_curr_stress.nc', 'Pair', 'Pair')




%% use stress right dx with just currents
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);

dx = (x_rho(1,2)-x_rho(1,1))/10;
dpdx = rho*ustar_av.ustrc.*ustar_av.ustrc;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_curr_stress.nc', 'Pair', 'Pair')









%% OTHER INPUT FILES



%% use stress right dx
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);

dx = (x_rho(1,2)-x_rho(1,1))/10;
dpdx = rho*ustar_av.ustrr.*ustar_av.ustrr;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_stress_dx.nc', 'Pair', 'Pair')
%% use stress right dx +30%
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);

dx = (x_rho(1,2)-x_rho(1,1))/10/1.3;
dpdx = rho*ustar_av.ustrr.*ustar_av.ustrr;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_stress_dx2.nc', 'Pair', 'Pair')
%%
ocean_time=ustar_av.dn-ustar_av.dn(1);

dx = 100;
dpdx = ustar_av.ustrr;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
%               dpdx = 0.02;
%               dpdx = 0.03*cos((ocean_time(it)*3600*24-tstart)*omega);
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs.nc', 'Pair', 'Pair')
%% use stress
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
dx = 100;
dpdx = rho*ustar_av.ustrr.*ustar_av.ustrr;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_stress.nc', 'Pair', 'Pair')

%% use stress/10
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
dx = 100;
fac=0.1;
dpdx = rho*ustar_av.ustrr.*ustar_av.ustrr*fac;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_stress_fac.nc', 'Pair', 'Pair')

%% use stress/12
ocean_time=ustar_av.dn-ustar_av.dn(1);
rho=1025;
dx = 100;
fac=1/12;
dpdx = rho*ustar_av.ustrr.*ustar_av.ustrr*fac;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_stress_fac2.nc', 'Pair', 'Pair')


%% use modeled ustrwm (m/s)
load('ustar_av')


ocean_time=ustar_av.dn-ustar_av.dn(1);

dx = 100;
dpdx = ustar_av.ustwm;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_wm.nc', 'Pair', 'Pair')
%% use modeled ustrc 
load('ustar_av')


ocean_time=ustar_av.dn-ustar_av.dn(1);

dx = 100;
dpdx = ustar_av.ustrc;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
%               dpdx = 0.02;
%               dpdx = 0.03*cos((ocean_time(it)*3600*24-tstart)*omega);
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_cur.nc', 'Pair', 'Pair')

%%
return
%% this is wrong! it is just the stress to velocity problem
clear

load('ustar_av')


ocean_time=ustar_av.dn-ustar_av.dn(1);

dx = 100;
fac=2/30;
dpdx = ustar_av.ustrr*fac;
dpdx(1)=dpdx(2);
for it=1:length(ocean_time)
    for j=1:7
        for i=1:8
           Pair(i,j,it)=1000+0.01*dpdx(it)*(i-1)*dx;
        end
    end
end

y_rho=linspace(-60,660,7);
x_rho=linspace(-41.6666666666667,541.666666666667,8);
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
create_roms_forcings_nogrid(x_rho',y_rho',ocean_time, 'Pair_coawst_velobs_adj.nc', 'Pair', 'Pair')

