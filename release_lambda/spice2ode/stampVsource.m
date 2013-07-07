function core=stampVsource(core, m, I, j, jb)
        if(j~=0),
            core.Y(m,j) = core.Y(m,j)+1;
            core.Y(j,I) = core.Y(j,I)+1;
        end
    
        if(jb~=0),
            core.Y(m,jb) = core.Y(m,jb)-1;
            core.Y(jb,I) = core.Y(jb,I)-1;
        end
end
