function drawDev(rrt_sample_t, tree, index)
    hold on
    figure(3)
    for i=2:size(tree,2)
        parent = tree(i).parent ;
        currentDev=tree(i).y(index)-getNominalValue(rrt_sample_t,  index, tree(i).t);
        parrentDev=tree(parent).y(index) - getNominalValue( rrt_sample_t, index, tree(parent).t);
        line( [tree(parent).t tree(i).t], [parrentDev currentDev] )
    end
    hold off
end