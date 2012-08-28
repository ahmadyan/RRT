function nearestNode = findNearestNode( tree, node )
    %this code is brute-force, needs optimization & vectorization.
    %todo: find nearest edge
    minimum_distance = 9909;
    nearestNode = 1 ; %root
    for i=1:size(tree,2)
        distance = sqrt ( ( node.x-tree(i).x)^2 + ( node.y-tree(i).y)^2 ) ;
        if (distance<minimum_distance),
            minimum_distance=distance ;
            nearestNode = i ;
        end
    end
end

