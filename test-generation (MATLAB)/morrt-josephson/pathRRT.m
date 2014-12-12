%% pathRRT
%%  - create a path from a start node to an end node
%%    using the RRT algorithm.
%%  - RRT = Rapidly-exploring Random Tree
%%  
%% Last Modified - 1/24/2006 - R. Beard

function pathRRT;

% create random world
Size = 100;
NumObstacles = 100;
world = createWorld(NumObstacles,[Size; Size],[0;0]);

% standard length of path segments
segmentLength = 5;

% randomly select start and end nodes
start_node = generateRandomNode(world);
end_node   = generateRandomNode(world);

% establish tree starting with the start node
tree = start_node;

% check to see if start_node connects directly to end_node
if ( (norm(start_node(1:2)-end_node(1:2))<segmentLength )...
    &(collision(start_node,end_node,world)==0) )
  path = [start_node; end_node];
else
  numPaths = 0;
  while numPaths<1,
      [tree,flag] = extendTree(tree,end_node,segmentLength,world);
      numPaths = numPaths + flag;
  end
end

% find path with minimum cost to end_node
path = findMinimumPath(tree,end_node);
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
    maxRadius = 5*maxRadius/NumObstacles/2;
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
chi      = 0;
cost     = 0;
node     = [pn, pe, chi, cost, 0];

% check collision with obstacle
while collision(node, node, world),
  pn       = (world.NEcorner(1)-world.SWcorner(1))*rand;
  pe       = (world.NEcorner(2)-world.SWcorner(2))*rand;
  chi      = 0;
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
    for sigma = 0:.2:1,
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
%% canEndConnectToTree
%%   check to see if the end node can connect to the tree
function flag = canEndConnectToTree(tree,end_node,minDist,world);
  flag = 0;
  % check only last node added to tree since others have been checked
  if ( (norm(tree(end,1:2)-end_node(1:2))<minDist)...
     & (collision(tree(end,1:2), end_node(1:2), world)==0) ),
    flag = 1;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extendTree
%%   extend tree by randomly selecting point and growing tree toward that
%%   point
function [new_tree,flag] = extendTree(tree,end_node,segmentLength,world);

  flag1 = 0;
  while flag1==0,
    % select a random point
    randomPoint = [...
        (world.NEcorner(1)-world.SWcorner(1))*rand,...
        (world.NEcorner(2)-world.SWcorner(2))*rand];

    % find leaf on node that is closest to randomPoint
    tmp = tree(:,1:2)-ones(size(tree,1),1)*randomPoint;
    [dist,idx] = min(diag(tmp*tmp'));
    cost     = tree(idx,4) + segmentLength;
    new_point = (randomPoint-tree(idx,1:2));
    new_point = tree(idx,1:2)+new_point/norm(new_point)*segmentLength;
    new_node = [new_point, 0, cost, idx];
    if collision(new_node, tree(idx,:), world)==0,
      new_tree = [tree; new_node];
      flag1=1;
    end
  end
  
  % check to see if new node connects directly to end_node
  if ( (norm(new_node(1:2)-end_node(1:2))<segmentLength )...
      &(collision(new_node,end_node,world)==0) )
    flag = 1;
    new_tree(end,3)=1;  % mark node as connecting to end.
  else
    flag = 0;
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% findMinimumPath
%%   find the lowest cost path to the end node
function path = findMinimumPath(tree,end_node);
    
    % find nodes that connect to end_node
    connectingNodes = [];
    for i=1:size(tree,1),
        if tree(i,3)==1,
            connectingNodes = [connectingNodes; tree(i,:)];
        end
    end

    % find minimum cost last node
    [tmp,idx] = min(connectingNodes(:,4));
    
    % construct lowest cost path
    path = [connectingNodes(idx,:); end_node];
    parent_node = connectingNodes(idx,5);
    while parent_node>1,
        parent_node = tree(parent_node,5);
        path = [tree(parent_node,:); path];
    end
    
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
      X = world.radius(i)*sin(th) + world.ce(i);
      Y = world.radius(i)*cos(th) + world.cn(i);
      fill(X,Y,'b');
  end
  
  X = path(:,2);
  Y = path(:,1);
  plot(X,Y,'r');


