clear
%% change magnitude of ripples
% !cp sed_floc_toy_ini.nc sed_floc_toy_ini_pre.nc
% !cp sed_floc_toy_ini.nc sed_floc_toy_ini_rip.nc
% !cp sed_floc_toy_ini.nc sed_floc_toy_ini_rip0.nc
%%
ncci=ncgeodataset('sed_floc_toy_ini.nc');
fno='sed_floc_toy_ini_rip.nc';
fno0='sed_floc_toy_ini_rip0.nc';

vars={'spherical'
    'Vtransform'
    'Vstretching'
    'theta_b'
    'theta_s'
    'Tcline'
    'hc'
    'Cs_r'
    'Cs_w'
    'sc_r'
    'sc_w'
    'ocean_time'
    'salt'
    'temp'
    'u'
    'ubar'
    'v'
    'vbar'
    'zeta'
    'mud_01'
    'mudfrac_01'
    'mudmass_01'
    'mud_02'
    'mudfrac_02'
    'mudmass_02'
    'mud_03'
    'mudfrac_03'
    'mudmass_03'
    'mud_04'
    'mudfrac_04'
    'mudmass_04'
    'mud_05'
    'mudfrac_05'
    'mudmass_05'
    'mud_06'
    'mudfrac_06'
    'mudmass_06'
    'mud_07'
    'mudfrac_07'
    'mudmass_07'
    'mud_08'
    'mudfrac_08'
    'mudmass_08'
    'mud_09'
    'mudfrac_09'
    'mudmass_09'
    'mud_10'
    'mudfrac_10'
    'mudmass_10'
    'mud_11'
    'mudfrac_11'
    'mudmass_11'
    'mud_12'
    'mudfrac_12'
    'mudmass_12'
    'mud_13'
    'mudfrac_13'
    'mudmass_13'
    'mud_14'
    'mudfrac_14'
    'mudmass_14'
    'mud_15'
    'mudfrac_15'
    'mudmass_15'
    'sand_01'
    'sandfrac_01'
    'sandmass_01'
    'bed_thickness'
    'bed_age'
    'bed_porosity'
    'bed_biodiff'
    'grain_diameter'
    'grain_density'
    'settling_vel'
    'erosion_stress'
    'ripple_length'
    'ripple_height'
    'dmix_offset'
    'dmix_slope'
    'dmix_time'
    'bed_tau_crit'};
   
nco=netcdf.open(fno,'NC_WRITE');
for ij=1:length(vars)
    dumid = netcdf.inqVarID(nco,vars{ij});
    if ij==[77]
        tempvar=ncci{vars{ij}}(:);
        tempvar=tempvar/10;
    else
        tempvar=ncci{vars{ij}}(:);
    end
    tempvar(isnan(tempvar))=nanmean(tempvar(:));
    netcdf.putVar(nco,dumid,tempvar);    
end
netcdf.close(nco);

nco=netcdf.open(fno0,'NC_WRITE');
for ij=1:length(vars)
    dumid = netcdf.inqVarID(nco,vars{ij});
    if ij==[77]
        tempvar=ncci{vars{ij}}(:);
        tempvar=tempvar*0;
    else
        tempvar=ncci{vars{ij}}(:);
    end
    tempvar(isnan(tempvar))=nanmean(tempvar(:));
    netcdf.putVar(nco,dumid,tempvar);    
end
netcdf.close(nco);




