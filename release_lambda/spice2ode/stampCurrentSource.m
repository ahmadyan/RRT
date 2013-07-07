function core=stampCurrentSource(core, i, j, x)
    core=stampJ(core, i, -x);
    core=stampJ(core, j, +x);
end