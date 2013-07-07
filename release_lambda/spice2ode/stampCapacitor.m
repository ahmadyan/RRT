function core=stampCapacitor(core, i)
    n1=core.Element(i).Node1;
    n2=core.Element(i).Node2;
    c= core.Element(i).Value;
    %Capacitor is open-circuited during dc analysis.
    if(core.mode==1),    %AC Analysis
        g= c*2*pi*j*core.freq;    
        %if(g>1e20) g=1e20;end
        core=stampConductance(core, n1, n2, g);  
        disp(core.Y)
        
    elseif(core.mode==2)    %Transient Analysis
        if(n1~=0), v1 = core.vNodeset(n1);
        else v1 = 0; end
        if(n2~=0), v2 = core.vNodeset(n2);
        else v2 = 0; end
            
        v=v1-v2;
	
        % capacitor companion model using trapezoidal approximation
	    % (Norton equivalent) consists of a current source in
	    % parallel with a resistor.  Trapezoidal is more accurate
	    % than backward euler but can cause oscillatory behavior
	    % if RC is small relative to the timestep.
		
        %g = core.dt/(2*c); %trapezoidal
        %curSourceValue = -v/g-current; 
		
        g = c/core.dt;     %BE
        i = -g*v;
        
        %   double voltdiff = volts[0] - volts[1];
	    %// we check compResistance because this might get called
	    %// before stamp(), which sets compResistance, causing
	    %// infinite current
	    %if (compResistance > 0)
		%current = voltdiff/compResistance + curSourceValue;
		%System.out.println("current=" + current + " voltdiff=" + voltdiff + " compresistance=" + compResistance + " curSourceValue=" + curSourceValue );
        
        %disp(sprintf('g=%f c=%f i=%f v=%f v1=%f v2=%f', g, c, i, v, v1, v2));
        core = stampConductance( core, n1, n2, g);
        core = stampCurrentSource(core, n1, n2, i );
    end
        
end
