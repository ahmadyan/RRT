function output = getNominalValue( rrt_sample_t, i, t  )
    x=floor( t/0.1e-4) +1;
    if(x>100)
        x=100;
    end
    output = rrt_sample_t(x, i);
end

