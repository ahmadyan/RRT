function config = generateConfig()
%UNTITLED3 handles program configuration
%   handles program configuration
global variation0
global variation1

%variables
config.iterations=3000;
config.learningSamples=300;
config.resolution=20;
config.sim_time=1e-4; ; %1
config.deltaT=1e-5; 
config.tmax=1e-3;
config.tbias=0.8;
config.numberOfClusters = 5 ;
config.segment = 0.01;
config.MAX=9909;
config.dim=15;
config.min=config.MAX*ones(1, config.dim);
config.max=-config.MAX*ones(1, config.dim);

variation0=0;
variation1=0;

disp('[config] running initial simulation to estimate min/max');
%{
t_end=config.tmax;
t=0;
steps=1000;
h=t_end/steps;
x=zeros(1,15);

for i=1:steps,
    x = ring(x, t, h);
	t=t+h;
    for j=1:config.dim,
            if (x(j)<config.min(j)),
                config.min(j) = x(j);
            end
            if (x(j)>config.max(j)),
                config.max(j) = x(j);
            end
    end
end

%add 10% variation to min/max
v=0.1;
for i=1:config.dim,
    config.min(i) = (1-sign(config.min(i))*v)*config.min(i);
    config.max(i) = (1+sign(config.max(i))*v)*config.max(i);
end
%}
config.min = [-0.6535   -0.5450   -0.7617   -0.7537   -0.7490   -0.7586   -0.9567   -3.3e-5  -3e-5    -0.0109   -0.0109   -0.0108   -0.0109   -0.0012   -0.0008];
config.max = [0.6499    0.5420    0.7522    0.7603    0.7573    0.7505    0.9568    0.5e-5    3e-5    0.0109    0.0109    0.0108    0.0108    0.0012    0.0008 ];
config.init = zeros(1, config.dim);

%tunnel-diode, depreched, moved to tunnel_diode.m
%x'=0.5*(y-(17.76*x-103.79*(x^2)+229.62*(x^3)-226.31*(x^4)+83.72*(x^5)))
%y'=0.2*(-x-1.5*y+1.2)
%-0.4:1.6;-0.4:1.6
%config.system = inline('[0.5*(y(2)-(17.76*y(1)-103.79*(y(1)^2)+229.62*(y(1)^3)-226.31*(y(1)^4)+83.72*(y(1)^5)));0.2*(-y(1)-1.5*y(2)+1.2)]','t','y');
%config.xmin=-4 ;
%config.xmax=+4 ;
%config.ymin=-4 ;
%config.ymax=+4 ;
%config.initX = -1  ;    % one of the three equilbrium point
%config.initY=  3  ;     % one of the three equilbrium point

config.generateRandomNodes=1;

config.sigma=0.05;
config.csigma=0;
config.alpha=0.5; %default value for alpha is 0.5
config.cluster=0; 
config.clusterSize=0;
config.maxcluster=10;

config.final_variance = 0.01 ;
config.init_variance = 0.5 ;
config.variance_cooling_rate = 1;

config.explore_step=100;

end

