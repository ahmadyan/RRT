function [core probe freq] = ac(config, core)
    core.mode=0;
    %initial DC Analysis to determine the operating points
    core=dc(config, core);
    core.mode=1; %switching to AC Analysis
    dcOperatingPoint = core.operatingPoint;
    disp('DC-Operating Condition:')
    dcOperatingPoint
    step=(core.fmax-core.fmin)/core.step;
    freq = core.fmin:step:core.fmax ;
    probe = zeros(size(freq,2),1);

    for i=1:size(freq,2),
        core.vNodeset = dcOperatingPoint(1:core.numNode) ;
        core.operatingPoint = dcOperatingPoint;
        disp(sprintf('AC Analysis for freq %f\n', freq(i)));
        core.freq=freq(i) ;
        core=stamp(config, core);
        core.operatingPoint = core.Y\core.J ;
        debug(core);
        probe(i) = core.operatingPoint(core.probe);
        
    end
    probe(1)=probe(2);
end