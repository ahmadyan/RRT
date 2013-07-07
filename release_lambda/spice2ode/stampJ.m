function core = stampJ(core, i, v);
    if(i~=0)
        core.J(i) = core.J(i) + v ;
    end
end