bed_thick = squeeze(ncread(url,'bed_thickness',[i j 1 1],[1 1 Inf Inf]));
bed_age = squeeze(ncread(url,'bed_age',[i j 1 1],[1 1 Inf Inf]));
ot=repmat(ocean_time,1,20)';
plot(mean(ot(2:10,:)-bed_age(2:10,:)))
hold on
plot(mean(ot(11:19,:)-bed_age(11:19,:)))