function config = generateConfig()
%UNTITLED3 handles program configuration
%   handles program configuration


%variables
config.iterations=300;
config.learningSamples=100;
config.resolution=20;
config.sim_time=10 ;
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

config.generateRandomNodes=0;
end

