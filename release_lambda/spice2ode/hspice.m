function [ flags ] = hspice( config )
%HSPICE hspice hook
%   Will connect to hspice and simulate the circuit.
%   TODO: bring back the hspice results into MATLAB

    if config.hspiceEnabled, 
        flags=system(sprintf('hspice %s', configfileName));
    else
        flags=-1;
    end

end

