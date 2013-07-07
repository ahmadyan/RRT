function core=stampConductance(core, n1, n2, g)
    %If neither side of the element is connected to ground
    %then subtract it from appropriate location in matrix.
    if (n1~=0) && (n2~=0),
        core.Y(n1,n2)= core.Y(n1,n2) - g;
        core.Y(n2,n1)= core.Y(n2,n1) - g;
    end
    
    %If node 1 is connected to ground, add element to diagonal of matrix.
    if (n1~=0),
        core.Y(n1,n1)= core.Y(n1,n1) + g;
    end
    %Ditto for node 2.
    if (n2~=0),
        core.Y(n2,n2)= core.Y(n2,n2) + g;
    end
end

