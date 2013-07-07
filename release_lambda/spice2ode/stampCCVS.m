function core=stampCCVS(core, m,k,kb,I,ID,r )
    if( k ~= 0 ),
        core.Y(m+1,k) = core.Y(m+1,k) +1 ;
        core.Y(k,I) = core.Y(k,I) +1 ;
    end
    if( kb ~= 0 ),
        core.Y(m+1,kb) = core.Y(m+1,kb) - 1 ;
        core.Y(kb,I) = core.Y(kb,I) - 1 ;
    end
    
    if( ID~=0),
        core.Y(m+1, ID) = core.Y(m+1,I) -r ;
    end
    
end