function core=stampVoltageSource(core, i)
    m=core.numNode+i;   %textbook writes this as m+1
    I=core.numNode+i;
    j=core.Vsource(i).Node1 ;
    jb=core.Vsource(i).Node2 ;
    v = core.Vsource(i).Value ;

    core=stampVsource(core, m, I, j, jb);
    if(core.mode==0)    %DC-Analysis
        if( core.Vsource(i).type==0)    %DC voltage source
            core=stampJ(core, m, v);    
        end        
    elseif (core.mode==1)   %AC-Analysis
        %if( core.Vsource(i).type==0)    %DC voltage source
        %    core=stampJ(core, m, v);    
        %else
        if(core.Vsource(i).type==1)
            a = core.Vsource(i).peak ;
            p = 0; %phase
            v = a +1j*p;
            core=stampJ(core, m, v); 
        end
    elseif (core.mode==2)   %Transient Analysis
        if( core.Vsource(i).type==0)    %DC voltage source
            core=stampJ(core, m, v);    
        elseif(core.Vsource(i).type==1)     %AC
            sim_time=core.step * core.dt;
            freqTimeZero=0;
            frequency=1;
            peak=core.Vsource(i).peak;
            offset=0;
            phaseShift=0;
            v = core.Vsource(i).peak ;
            %u = v * exp (-d * t * f) * sin (2 * pi * f * t + rad (p));
            %ignoring the ac sources during the transient analysis for now.
  
        elseif(core.Vsource(i).type==2) %VSIN
            sim_time=core.step * core.dt;
            freqTimeZero=0;
            frequency=core.Vsource(i).freq;
            peak=core.Vsource(i).peak;
            offset=core.Vsource(i).offset;
            phaseShift=0;
            w = 2*pi*(sim_time-freqTimeZero)*frequency + phaseShift;
            v = sin(w)*peak+offset;
            core=stampJ(core, m, v);    
        end
    end
end

% 	case WF_DC: return maxVoltage+bias;
% 	case WF_AC: return Math.sin(w)*maxVoltage+bias;
% 	case WF_SQUARE: bias+((w % (2*pi) > (2*pi*dutyCycle)) ? -maxVoltage : maxVoltage);
% 	case WF_TRIANGLE: bias+triangleFunc(w % (2*pi))*maxVoltage;
% 	case WF_SAWTOOTH: bias+(w % (2*pi))*(maxVoltage/pi)-maxVoltage;
% 	case WF_PULSE:((w % (2*pi)) < 1) ? maxVoltage+bias : bias;
	
    