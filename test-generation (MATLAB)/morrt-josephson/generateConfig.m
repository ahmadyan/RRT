function config = generateConfig()
%UNTITLED3 handles program configuration
%   handles program configuration


%variables
config.iterations=3000;
config.learningSamples=300;
config.resolution=20;
config.sim_time=10 ; %1
config.deltaT=0.1  ;
config.numberOfClusters = 5 ;
config.segment = 0.01;
config.MAX=9909;

%tunnel-diode, depreched, moved to tunnel_diode.m
%x'=0.5*(y-(17.76*x-103.79*(x^2)+229.62*(x^3)-226.31*(x^4)+83.72*(x^5)))
%y'=0.2*(-x-1.5*y+1.2)
%-0.4:1.6;-0.4:1.6
config.system = inline('[0.5*(y(2)-(17.76*y(1)-103.79*(y(1)^2)+229.62*(y(1)^3)-226.31*(y(1)^4)+83.72*(y(1)^5)));0.2*(-y(1)-1.5*y(2)+1.2)]','t','y');
config.xmin=-4 ;
config.xmax=+4 ;
config.ymin=-4 ;
config.ymax=+4 ;
config.initX = -1  ;    % one of the three equilbrium point
config.initY=  3  ;     % one of the three equilbrium point


config.generateRandomNodes=1;

config.sigma=0.01;
config.csigma=0;
config.alpha=0.8; %default value for alpha is 0.5
config.cluster=0; 
config.clusterSize=0;
config.maxcluster=10;

config.final_variance = 0.01 ;
config.init_variance = 0.2 ;
config.variance_cooling_rate = 1;
end

