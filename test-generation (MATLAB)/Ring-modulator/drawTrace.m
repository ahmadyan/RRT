function drawTrace(tree, index)
    hold on
    figure(2)
    for i=2:size(tree,2)
        parent = tree(i).parent ;
        line( [tree(parent).t tree(i).t], [tree(parent).y(index) tree(i).y(index)] )
    end
    hold off
end