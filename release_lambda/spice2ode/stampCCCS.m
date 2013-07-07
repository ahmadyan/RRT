function core=stampCCCS(core, k, kb, I, ID, a, m)
  %stamping the dependant voltage source
    if( k ~= 0 ),
        core.Y(k,ID) = core.Y(k,ID) +a ;
    end
    
    if( kb ~= 0 ),
        core.Y(kb,ID) = core.Y(kb,ID) - a ;
    end
    
    %I1=a*I2:
    
    core.Y(m+1,ID)=core.Y(m+1,ID)-a;
    core.Y(m+1,I)=core.Y(m+1,I)+1;
end