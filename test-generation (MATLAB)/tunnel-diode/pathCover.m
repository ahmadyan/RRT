%% pathCover
%%  - create a path that covers a region with obstacles
%%  
%% Last Modified - 1/24/2006 - R. Beard

function pathCover;

% create random world
Size = 100;
NumObstacles = 10;
world = createWorld(NumObstacles,[Size; Size],[0;0]);

% model of world return map
mapSize = 10;  % this is a critical parameter!
return_map = 50*ones(mapSize,mapSize)+ rand(mapSize,mapSize);
%return_map = 50*ones(mapSize,mapSize);
%plotReturnMap(return_map), pause

% minimum turn radius 
turnRadius = 5;

% randomly generate initial node
node0=generateRandomNode(world);
path = node0;

for i=1:200
  % generate tree from current path node
  tree = generateTree(path(end,:),world,turnRadius,return_map);

  % find path with highest return and return next node
  node = highestReturnPath(tree);
  
  % modify return map to indicate visited location
  return_map = updateReturnMap(path(end,:), return_map, world);
%  plotReturnMap(return_map), pause  
  
  % add node to the path
  path = [path; node];
  
end

plotWorld(world,path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% createWorld
%%  - create random world with obstacles
%%  the first element is the north coordinate
%%  the second element is the south coordinate
function world = createWorld(NumObstacles, NEcorner, SWcorner);

  % check to make sure that the region is nonempty
  if (NEcorner(1) <= SWcorner(1)) | (NEcorner(2) <= SWcorner(2)),
      disp('Not valid corner specifications!')
      world=[];
      
  % create world data structure
  else
    world.NumObstacles = NumObstacles;
    world.NEcorner = NEcorner;
    world.SWcorner = SWcorner;
                          
    % create NumObstacles 
    maxRadius = min(NEcorner(1)- SWcorner(1), NEcorner(2)-SWcorner(2));
    maxRadius = maxRadius/NumObstacles/2;
    for i=1:NumObstacles,
        % randomly pick radius
        world.radius(i) = maxRadius*rand;
        % randomly pick center of obstacles
        cn = SWcorner(1) + world.radius(i)...
            + (NEcorner(1)-SWcorner(1)-2*world.radius(i))*rand;
        ce = SWcorner(2) + world.radius(i)...
            + (NEcorner(2)-SWcorner(2)-2*world.radius(i))*rand;
        world.cn(i) = cn;
        world.ce(i) = ce;
    end
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generateRandomNode
%%   create a random node (initialize)
function node=generateRandomNode(world);

% randomly pick configuration
pn       = (world.NEcorner(1)-world.SWcorner(1))*rand;
pe       = (world.NEcorner(2)-world.SWcorner(2))*rand;
chi      = 2*pi*rand;
cost     = 0;
node     = [pn, pe, chi, cost, 0];

% check collision with obstacle
while collision(node, node, world),
  pn       = (world.NEcorner(1)-world.SWcorner(1))*rand;
  pe       = (world.NEcorner(2)-world.SWcorner(2))*rand;
  chi      = 2*pi*rand;
  cost     = 0;
  node     = [pn, pe, chi, cost, 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% collision
%%   check to see if a node is in collsion with obstacles
function collision_flag = collision(node, parent, world);

collision_flag = 0;

if ((node(1)>world.NEcorner(1))...
    | (node(1)<world.SWcorner(1))...
    | (node(2)>world.NEcorner(2))...
    | (node(2)<world.SWcorner(2)))
  collision_flag = 1;
else
    for sigma = 0:.5:1,
    p = sigma*node(1:2) + (1-sigma)*parent(1:2);
      % check each obstacle
      for i=1:world.NumObstacles,
        if (norm([p(1);p(2)]-[world.cn(i); world.ce(i)])<=1.5*world.radius(i)),
            collision_flag = 1;
            break;
        end
      end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generateTree
%%   Generate the search Tree (three levels
function tree = generateTree(node0,world,turnRadius,return_map);

  nodes_expanded = 0;
  tree = node0;
  while nodes_expanded<3^0+3^1+3^2,  % tree depth is three
    nodes_expanded = nodes_expanded+1;
    new_nodes = expandNode(tree,nodes_expanded,turnRadius,world,return_map);
    tree = [tree; new_nodes];
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% expandNode
%%   expand a node to produce three next step nodes 
function new_nodes = expandNode(tree,parent,turnRadius,world,return_map);

  parent_node = tree(parent,:);
  new_nodes = [];
  % turn right pi/2
  p = [parent_node(1); parent_node(2)]...
      + turnRadius*[-sin(parent_node(3))+cos(parent_node(3));...
                     cos(parent_node(3))+sin(parent_node(3))];
  chi  = parent_node(3)+pi/2;
  cost = parent_node(4) + findReturn(p(1),p(2),return_map,world);
  node = [p(1), p(2), chi, cost, parent];
  if collision(node, parent_node, world)==0,
      new_nodes = [new_nodes; node];
  end
      
  
  % go straight
  % turn right pi/2
  p = [parent_node(1); parent_node(2)]...
      + pi/2*turnRadius*[cos(parent_node(3)); sin(parent_node(3))];
  chi  = parent_node(3);
  cost = parent_node(4) + findReturn(p(1),p(2),return_map,world);
  node = [p(1), p(2), chi, cost, parent];
  if collision(node, parent_node, world)==0,
      new_nodes = [new_nodes; node];
  end

  % turn left pi/2
  p = [parent_node(1); parent_node(2)]...
      + turnRadius*[sin(parent_node(3))+cos(parent_node(3));...
                    -cos(parent_node(3))+sin(parent_node(3))];
  chi  = parent_node(3)-pi/2;
  cost = parent_node(4) + findReturn(p(1),p(2),return_map,world);
  node = [p(1), p(2), chi, cost, parent];
  if collision(node, parent_node, world)==0,
      new_nodes = [new_nodes; node];
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% highestReturnPath
%%   find path with highest return and return next node
function node = highestReturnPath(tree);

first_leaf_node = 3^0+3^1+3^2+1; % tree depth is tree
[foo,idx_tmp] = max(tree(first_leaf_node:end,4));
idx = first_leaf_node-1+idx_tmp;

idx_parent = tree(idx,5);
while idx_parent > 1,
    idx = idx_parent;
    idx_parent = tree(idx,5);
end
    
node = tree(idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% updateReturnMap
%%   update the return map to indicate where MAV has been
function new_return_map = updateReturnMap(node,return_map,world);

  new_return_map = return_map;
  pn = node(1);
  pe = node(2);
  [pn_max,pe_max] = size(return_map);
  fn = pn_max*(pn/(world.NEcorner(1)-world.SWcorner(1)));
  fn = min(pn_max,round(fn));
  fn = max(1,fn);
  fe = pe_max*(pe/(world.NEcorner(2)-world.SWcorner(2)));
  fe = min(pe_max,round(fe));
  fe = max(1,fe);
  
  new_return_map(fn,fe) = return_map(fn,fe) - 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% findReturn
%%   compute the return value at a particular location
function return_value = findReturn(pn,pe,return_map,world);

  [pn_max,pe_max] = size(return_map);
  fn = pn_max*(pn/(world.NEcorner(1)-world.SWcorner(1)));
  fn = min(pn_max,round(fn));
  fn = max(1,fn);
  fe = pe_max*(pe/(world.NEcorner(2)-world.SWcorner(2)));
  fe = min(pe_max,round(fe));
  fe = max(1,fe);
  return_value = return_map(fn,fe);
  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotWorld
%%   plot obstacles and path
function plotWorld(world,path)
  % the first element is the north coordinate
  % the second element is the south coordinate

  N = 10;
  th = [0:2*pi/N:2*pi];
  figure(1), clf
  axis([world.SWcorner(2),world.NEcorner(2),...
      world.SWcorner(1), world.NEcorner(1)]);
  hold on
  
  for i=1:world.NumObstacles,
      X = world.radius(i)*cos(th) + world.ce(i);
      Y = world.radius(i)*sin(th) + world.cn(i);
      fill(X,Y,'b');
  end
  
  X = path(:,2);
  Y = path(:,1);
  plot(X,Y,'r');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotReturnMap
%%   plot the return map
function plotReturnMap(return_map);

figure(2), clf 
mesh(return_map)  
