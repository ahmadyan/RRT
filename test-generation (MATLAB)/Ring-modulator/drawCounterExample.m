function drawCounterExample( tree )
%drawTree Draws a Tree
%   Detailed explanation goes here
hold on
figure(5)
node = tree( size(tree,2) );
while( node.parent ~= -1 ),   
     line([node.x tree(node.parent).x ], [ node.y tree(node.parent).y] )
      node = tree( node.parent ) ;
end


hold off

end

