function config=configuration(fileName)
    %config.fileName = 'circuit-rc1.txt' ;
    %config.fileName = 'two_section_rc.txt'; %for ac-testign
    %config.fileName = 'amp_tran1.txt' ;
    %config.fileName  = '741_amplifier.txt' ;
    config.fileName = fileName ;
    %config.fileName = 'diode.txt' ;
    %config.fileName = 'bridge-t.txt'
    config.verbose=0;
    config.symbolic=0;
    config.hspiceEnabled=0;
    
    config.DCErrorThreshold = 10^(-7);
    
end