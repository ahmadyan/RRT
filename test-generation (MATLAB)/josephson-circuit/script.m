%AUTHOR: Adel
%DATE: Nov. 25 (BLKFRD)


config = generateConfig();
cluster=0; clusterSize=0;
disp('doing data-mining stuff');
if config.generateRandomNodes==0,
    [cluster clusterSize] = datamining(config);
end
disp('Generating node sequence');
nodes = generateNewNode(config, cluster, clusterSize);
disp('All-done, now generating RRT')
tree = generateTree(config, nodes);
disp('finished, drawing the tree');
drawTree(tree);


%do_data_mining();
%create_tree();
%expand_tree();
%generate_report();

%hold on
%should be replaced with generatingRootCode
%root = generateNewNode(config);
%scatter(root)
%for i=1:(config.iterations),
%    node = generateNewNode(config);
    %scatter(node)
    %q = findNearestNode( tree, node ) ;
 %   plot(root, node )
    %addEdge( tree, q, node );
%end

tree2 = generateTree2( config, tree, 0.75, 0.5, 50);
tree3 = generateTree2( config, tree2, 0, 0.5, 50) ;
tree4 = generateTree2( config, tree3, -0.5, 0, 50) ;
tree5 = generateTree2( config, tree4, -0.75, 0.25, 50) ;
tree6 = generateTree2( config, tree5, 0, 0, 50) ;
tree6 = generateTree2( config, tree6, 0, 0, 100) ;
drawTree(tree6)