function drawTree( tree )
%drawTree Draws a Tree
%   Detailed explanation goes here
hold on
figure(4)
for i=2:size(tree,2)
    parent = tree(i).parent ;
    %line([tree(i).x tree(i).x ], [ tree(i).y tree(i).y], 'Marker','.','LineStyle','-')
    line([tree(i).x tree(parent).x ], [ tree(i).y tree(parent).y] )
    %line([tree(i).x tree(i).x+0.01], [ tree(i).y tree(i).y+0.01] )
    tree(i).x ;
    
end
hold off

end

