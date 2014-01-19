function drawTree( tree, index1, index2 )
%drawTree Draws a Tree
%   Detailed explanation goes here
hold on
figure(4)
for i=2:size(tree,2)
    parent = tree(i).parent ;
    %line([tree(i).x tree(i).x ], [ tree(i).y tree(i).y], 'Marker','.','LineStyle','-')
    %line([tree(i).y(index1) tree(parent).y(index1) ], [ tree(i).y(index2) tree(parent).y(index2)] )
    line([tree(i).y(index1) tree(i).y(index1)+0.01], [ tree(i).y(index2) tree(i).y(index2)+0.01] )
end
hold off

end

