function core = stampVCVS(core, m, I, k, kb, j, jb, u )        
    if( j ~= 0 ),
        core.Y(m+1,j) = core.Y(m+1,j) - u ;
    end
    
    if( jb ~= 0 ),
        core.Y(m+1,jb) = core.Y(m+1,jb) + u ;
    end
    
    if( k ~= 0 ),    
        core.Y(m+1,k) = core.Y(m+1,k) +1 ;
        core.Y(k,I) = core.Y(k,I) +1 ;
    end
    
    if( kb ~= 0 ),
        core.Y(m+1,kb) = core.Y(m+1,kb) - 1 ;
        core.Y(kb,I) = core.Y(kb,I) - 1 ;
    end

end