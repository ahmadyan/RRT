function [core probe probe2 probe3] = tr(config, core)
    core.mode=0;
    core=dc(config, core);
    Y=core.Y;
    J=core.J;
    core.mode=2;    
    disp('-------------DC-COMPLETE-------------');
    totalSteps = floor( (core.tf - core.t0) / core.dt );
    probe=zeros(totalSteps, 1);
    probe2=zeros(totalSteps, 1);
    probe3=zeros(totalSteps, 1);
    for step=1:totalSteps,
        core.step=step;
        core=dc(config, core);
%         disp('transient');
%         disp('Y');
%         core.Y
%         disp('J');
%         core.J
%         disp('X');
%         core.operatingPoint
        probe(step)=core.operatingPoint(core.probe);
        probe2(step)=core.operatingPoint(core.probe2);
        probe3(step)=core.operatingPoint(core.probe3);
    end
end