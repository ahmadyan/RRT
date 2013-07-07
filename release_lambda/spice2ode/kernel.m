function [config core] = kernel(fileName, mode, probe1, probe2, probe3)
tic
%clear;clc;
disp(sprintf('[info] spice started!'));
config = configuration(fileName);
config.fileName=fileName;
disp(sprintf('[info] parsing input netlist: %s', config.fileName)); 

hspice(config);
core = init(mode, probe1, probe2, probe3);
core=parser(config, core);
switch core.mode
    case 0,
        core=dc( config, core );
    case 1,
        [core probe freq]=ac( config, core );
    case 2,
        [core probe, probe2, probe3]=tr(config,core);
end
%generating report
if(core.mode==0)
    core.operatingPoint
elseif(core.mode==1)
    reportGenerator(config, core, freq, probe);
    probe
else
    figure(1); plot(probe);
    figure(2); plot(probe2);
    figure(3); plot(probe3);
end
disp(sprintf('[info] Elapsed time = %g seconds.\n',toc));
end