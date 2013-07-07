function core = StampVCCS( core, g, k, kb, j, jb )
    if ((j~=0)&&(k~=0)),
        core.Y(k,j) = core.Y(k,j)+g;
    end
    if ((j~=0)&&(kb~=0)),
        core.Y(kb,j) = core.Y(kb,j)-g;
    end
    if ((jb~=0)&&(k~=0)),
        core.Y(k,jb) = core.Y(k,jb)-g;
    end
    if ((jb~=0)&&(kb~=0)),
        core.Y(kb,jb) = core.Y(kb,jb)+g;
    end
end

