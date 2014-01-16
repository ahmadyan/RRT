tree2 = generateTree2( config, tree2, 1, 0.1, 100);
tree2 = generateTree2( config, tree, 0.75, 0.5, 100);
tree2 = generateTree2( config, tree2, 0.5, 0.75, 100);
tree2 = generateTree2( config, tree2, 0, 1, 100);
tree2 = generateTree2( config, tree2, -0.25, 0.75, 100);
tree2 = generateTree2( config, tree2, -0.5, 0.5, 200);
disp('x');
tree2 = generateTree2( config, tree2, -0.25, -0.25, 200);
tree2 = generateTree2( config, tree2, 0.25, -0.25, 200);
disp('x');
tree2 = generateTree2( config, tree2, 0.25, -0.25, 200);
tree2 = generateTree2( config, tree2, -0.25, 0.25, 200);
disp('x');
tree2 = generateTree2( config, tree2, 0, 0, 800);

drawTree(tree2);