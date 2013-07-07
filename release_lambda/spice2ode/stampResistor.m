function core=stampResistor(core, i)
    n1=core.Element(i).Node1;
    n2=core.Element(i).Node2;
    g= 1/core.Element(i).Value;
	core=stampConductance(core, n1, n2, g);                
end

