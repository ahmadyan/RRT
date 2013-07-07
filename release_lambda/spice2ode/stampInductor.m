function core=stampInductor(core, i)
    if(core.ac==1)
        n1=core.Element(i).Node1;
        n2=core.Element(i).Node2;
        l= core.Element(i).Value;
        g= 1/(l*2*pi*j*core.freq);
        core=stampConductance(core, n1, n2, g);  
    end
        
end
