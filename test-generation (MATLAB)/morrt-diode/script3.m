%name: script3
%date: jan-28th-2012
%objective: 1. generate a classic RRT algorithm
%           2. observe the growth of RRT as algorithm progresses

%step 1: generate a classic RRT 
config = generateConfig();
size=100;
cluster=0; clusterSize=0; maxcluster=10;
nodes = generateNewNode(config, cluster, clusterSize);
tree = generateTree(config, nodes);

sample = zeros(2, size);
for i=1:size,
    sample(1, i) = tree(i).x;
    sample(2, i) = tree(i).y; 
end
drawTree(tree);

figure(2)
label = vbgm(sample, maxcluster); 
spread(sample, label)