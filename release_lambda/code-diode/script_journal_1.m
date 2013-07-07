%{
clc; clear;
config = generateConfig();
config.alpha=0.01;
script
save rrt_a_0_0.mat



clc; clear;
config = generateConfig();
config.alpha=0.1;
script
save rrt_a_0_1.mat


clc; clear;
config = generateConfig();
config.alpha=0.2;
script
save rrt_a_0_2.mat

clc; clear;
config = generateConfig();
config.alpha=0.3;
script
save rrt_a_0_3.mat


clc; clear;
config = generateConfig();
config.alpha=0.4;
script
save rrt_a_0_4.mat



clc; clear;
config = generateConfig();
config.alpha=0.5;
script
save rrt_a_0_5.mat

clc; clear;
config = generateConfig();
config.alpha=0.6;
script
save rrt_a_0_6.mat

clc; clear;
config = generateConfig();
config.alpha=0.7;
script
save rrt_a_0_7.mat

clc; clear;
config = generateConfig();
config.alpha=0.8;
script
save rrt_a_0_8.mat

clc; clear;
config = generateConfig();
config.alpha=0.9;
script
save rrt_a_0_9.mat
%}
clc; clear;
config = generateConfig();
config.alpha=0.99;
script
save rrt_a_1_0.mat