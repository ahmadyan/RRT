function config = generateConfig()
%UNTITLED3 handles program configuration
%   handles program configuration


%variables
config.iterations=3000;
config.learningSamples=100;
config.resolution=20;
config.sim_time=0.7 ;
config.deltaT=0.7;
config.numberOfClusters = 5 ;
config.segment = 0.01;
config.MAX=9909;

%tunnel-diode, depreched, moved to tunnel_diode.m
%x'=0.5*(y-(17.76*x-103.79*(x^2)+229.62*(x^3)-226.31*(x^4)+83.72*(x^5)))
%y'=0.2*(-x-1.5*y+1.2)
%-0.4:1.6;-0.4:1.6
config.system_forward = inline('[0.5*(y(2)-(17.76*y(1)-103.79*(y(1)^2)+229.62*(y(1)^3)-226.31*(y(1)^4)+83.72*(y(1)^5)));0.2*(-y(1)-1.5*y(2)+1.2)]','t','y');

config.system_backward = inline('[0.5*(y(2)-(17.76*y(1)-103.79*(y(1)^2)+229.62*(y(1)^3)-226.31*(y(1)^4)+83.72*(y(1)^5)));0.2*(-y(1)-1.5*y(2)+1.2)]','-t','y');

config.xmin=-0.4 ;
config.xmax=+1.6 ;
config.ymin=-0.4 ;
config.ymax=+1.6 ;
config.initX=0.285 ;    % one of the three equilbrium point
config.initY=0.61 ;     % one of the three equilbrium point

config.generateRandomNodes=0;
end


