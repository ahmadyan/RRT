function nearestNode = findNearestNode( tree, node , config)
    %this code is brute-force, needs optimization & vectorization.
    %todo: find nearest edge
    minimum_distance = 9909;
    nearestNode = 1 ; %root
    
    for i=1:config.i-1,
        distance=0;
        for j=1:size(node.y, 2),
            distance = distance + ( (node.y(j)-tree(i).y(j))/(config.max(j)-config.min(j)) )^2;
        end
        distance = sqrt(distance);
        
        %tdistance = abs( node.t - tree(i).t ) / config.tmax ;
        %distance = ((1-config.tbias)*distance + config.tbias*tdistance) ;
        if (distance<minimum_distance),
            minimum_distance=distance ;
            nearestNode = i ;
        end
    end
end

