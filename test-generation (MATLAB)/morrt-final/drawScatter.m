function drawScatter( tree, index1, index2 )
%drawTree Draws a Tree
%   Detailed explanation goes here
hold on
figure(4)
x = zeros(1, size(tree, 2));
y = zeros(1, size(tree, 2));
for i=1:size(tree,2)
    %parent = tree(i).parent ;
    x(i) = tree(i).y(index1); 
    y(i) = tree(i).y(index2);
    
end
scatter( x, y);
hold off

end

