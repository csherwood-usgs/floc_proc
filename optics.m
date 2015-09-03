% optics.m
% instrument response for our fractal dimension and floc diameters
load c_coeffs
load c_coeffs_sand

c_ac9_nf = interp1(nf',c_ac9',fnf)';
c_ac9_i = interp1(D*1e-6',c_ac9_nf,fdiam')';
r_ac9=sum(C.*repmat(c_ac9_i,nz,1),2);

c_mass_nf = interp1(nf',c_mass',fnf)';
c_mass_i = interp1(D*1e-6',c_mass_nf,fdiam')';
r_mass=sum(C.*repmat(c_mass_i,nz,1),2);

c_cstar_nf = interp1(nf',c_cstar',fnf)';
c_cstar_i = interp1(D*1e-6',c_cstar_nf,fdiam')';
r_cstar=sum(C.*repmat(c_cstar_i,nz,1),2);

c_lisst_nf = interp1(nf',c_lisst',fnf)';
c_lisst_i = interp1(D*1e-6',c_lisst_nf,fdiam')';
r_lisst=sum(C.*repmat(c_lisst_i,nz,1),2);

% optical instrument response for sand diameters
Sc_ac9_i = interp1(D*1e-6',c_ac9_sand,Dm')';
Sr_ac9=sum(SC.*repmat(Sc_ac9_i,nz,1),2);

Sc_mass_i = interp1(D*1e-6',c_mass_sand,Dm')';
Sr_mass=sum(SC.*repmat(Sc_mass_i,nz,1),2);

Sc_cstar_i = interp1(D*1e-6',c_cstar_sand,Dm')';
Sr_cstar=sum(SC.*repmat(Sc_cstar_i,nz,1),2);

Sc_lisst_i = interp1(D*1e-6',c_lisst_sand,Dm')';
Sr_lisst=sum(SC.*repmat(Sc_lisst_i,nz,1),2);
r_ac9 = r_ac9+Sr_ac9;
r_mass = r_mass+Sr_mass;
r_cstar = r_cstar+Sr_cstar;
r_lisst = r_lisst+Sr_lisst;