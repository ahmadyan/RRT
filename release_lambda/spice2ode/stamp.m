function core = stamp(config, core)
    % Preallocate all of the cell arrays
    core.Y=zeros(core.numNode + core.numV + core.numE + core.numF + core.numH , core.numNode + core.numV +  core.numE + core.numF + core.numH);
    %core.Y=sparse(core.numNode + core.numV + core.numE + core.numF + core.numH , core.numNode + core.numV +  core.numE + core.numF + core.numH);
    core.V=cell(core.numNode+ core.numV + core.numE + core.numF + core.numH,1);
    core.J=zeros(core.numNode+ core.numV + core.numE + core.numF + core.numH,1);
    %core.J=sparse(core.numNode+ core.numV + core.numE + core.numF + core.numH,1);

    %stamping elements
    for i=1:core.numElem,
       switch(core.Element(i).Name(1)),
            case {'R' , 'r'},
                core=stampResistor(core, i);
            case {'L', 'l'},
            	core=stampInductor(core, i);
            case {'C', 'c'},
                core=stampCapacitor(core, i);
            case {'D', 'd'},
                core=stampDiode(core, i);
            case {'Q', 'q'},
            	core=stampBJT(core, i);
       end
    end

    %stamping sources (current & voltages)
    for i=1:core.numI,
         core=stampCurrentSource(core, core.Isource(i).Node1, core.Isource(i).Node2, core.Isource(i).Value);
    end
    
    %Voltage-Controlled Voltage Sources
    for i=1:core.numV,
        core=stampVoltageSource(core, i);
        
    end

    for i=1:core.numE,
        v = 0 ;
        m = core.numNode+ core.numV + i -1 ;
        I = core.numNode+ core.numV + i ;
        k = core.Esource(i).Node1 ;
        kb= core.Esource(i).Node2 ;
        j = core.Esource(i).Node3 ;
        jb= core.Esource(i).Node4 ;  
        u = core.Esource(i).Value ;    
        core=stampVCVS(core, m, I, k, kb, j, jb, u);
        core=stampJ(core, m+1, v);
    end
    
    %Linear Current-Controlled Current Sources
    for i=1:core.numF,
        x = core.numNode+core.numV+core.numE+i ;
        v = 0 ;  % Vsource( Fsource(i).DependantID ).Value ??
        
        m = core.numNode+core.numV + core.numE + i -1 ;
        I= core.numNode+core.numV + core.numE +i ;
        ID = core.numNode + core.Fsource(i).DependantID ;
        k=core.Fsource(i).Node1 ;
        kb=core.Fsource(i).Node2 ;
        j=core.Fsource(i).Node3 ;
        jb=core.Fsource(i).Node4 ;  
        a=core.Fsource(i).Value ;
    
        core=stampCCCS(core, k, kb, I, ID, a, m); 
        core=stampJ(core, x, v);
    end
 
     
    %Linear Voltage-Controlled Current Sources
    for i=1:core.numG,
        k=core.Gsource(i).Node1 ;
        kb=core.Gsource(i).Node2 ;
        j=core.Gsource(i).Node3 ;
        jb=core.Gsource(i).Node4 ;  
        g=core.Gsource(i).Value ;
        core=StampVCCS( core, g, k, kb, j, jb );
    end


    %Linear Current-Controlled Voltage Sources:
    for i=1:core.numH,
        x = core.numNode+core.numV+core.numE+core.numF+ i ;
        v = 0 ;
        m = core.numNode+core.numV + core.numE + core.numF + i -1 ;
        I= core.numNode+core.numV + core.numE + core.numF + i  ;
        ID = core.numNode + core.Hsource(i).DependantID ;
        k=core.Hsource(i).Node1 ;
        kb=core.Hsource(i).Node2 ;
        j=core.Hsource(i).Node3 ;
        jb=core.Hsource(i).Node4 ;  
        r=core.Hsource(i).Value ;
        
        core=stampJ(core, x, v);
        core=stampCCVS(core, m,k,kb,I,ID,r );
    end
    
end    