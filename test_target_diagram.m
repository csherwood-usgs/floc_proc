% script to test target diagram
figure(1);clf;
data = cos(0:2*pi/20:4*pi);
model = 1+data; % bias high
target_diagram(model,data,1)
model = 2*data; % more variability
target_diagram(model,data,1)
model = rand(size(data))+cos(pi/4+(0:2*pi/20:4*pi)); % noise and out of phase
target_diagram(model,data,1)