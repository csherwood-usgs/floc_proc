function [f, chi] = gfp(a)
% [f,chi]=gpf(a) - Get floc functions from Thorne plot for a (microns)
% The floc fractal dimention is 2.0 and g_o = h_o = 1.01, whatever that
% means.

% This is CRS version; slower than ALA version, but better match
% at ends. This is an attempt to match the paper version faxed from Thorne
% and has been replaced by f_chi_func.m
% Compare with get_floc_param.m
% csherwood@usgs.gov 22 Feb 2014
af = [.2  2 200 1000 2000];
fe = [-6 -4 -1.785 -1.785 -1.785];
ac = [.1 .7 2 50 1000 2000];
chie = [-5 -3.2 -3 -5 -2.77 -2.4];
f = 10.^interp1(log10(af),fe,log10(a),'linear');
chi = 10.^interp1(log10(ac),chie,log10(a),'linear');