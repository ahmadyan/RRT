function core = stampY( core, i, j, v )
    if(i~=0&j~=0)
        core.Y(i,j) = core.Y(i,j) + v;
    end
end

