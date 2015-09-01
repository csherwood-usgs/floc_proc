% optics_response - Response to various grain sizes and density
% area of flocs
d= fdiam(1:NST);
r= rhos(1:NST);
Afe =  (3./4.)*( .1 ./r)./(0.5*d);

ds=fdiam(NST+NNN+1:NST+NNN+NND);
rs = rhos(NST+NNN+1:NST+NNN+NND);
Afs = (3./4.)*( .1 ./rs)./(0.5*ds);