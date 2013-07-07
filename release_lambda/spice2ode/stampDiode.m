function core=stampDiode(core, i)
    n1=core.Element(i).Node1;
    n2=core.Element(i).Node2;
    if(n1~=0)
        v1 = core.vNodeset(n1);
    else
        v1 = 0;
    end
    
    if(n2~=0)
        v2 = core.vNodeset(n2);
    else
        v2 = 0;
    end
            
	%diode equations:
    Is = 10^(-15);  % saturation current
    Vt = 0.02585126075 ;    % 1/40  % thermal voltage
    
	vd=v1-v2;
	id = Is * ( exp(vd/Vt)-1 ); %standard diode equation
	geq = id/Vt ;       % g_eq = dId/dVd = IS/Vt * exp(Vd/vt) = i(v_n)/V_t         
    Ieq = id - geq*vd ;
    
    core=stampConductance(core, n1, n2, geq);
    core=stampCurrentSource(core, n1, n2, Ieq);
    core.isNonlinear=1;        
end

