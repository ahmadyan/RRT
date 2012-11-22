% Rapidly-Exploring Random Forests
% Adel Ahmadyan, Nov-2012
clear
clc
%setup
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
xmin=-0.3 ;
xmax=+1.3 ;
ymin=-0.3 ;
ymax=+1.3 ;
xinit= 0.285;
yinit= 0.61;
deltaT=0.1;
  
num_roots = 23 ;
num_nodes = 5000;
%generate random points
for i=1:num_nodes,
    % randomly pick configuration
    node.x = (xmax-xmin)*rand + xmin;
    node.y = (ymax-ymin)*rand + ymin;
    node.parent  = -1 ;
    node.direction = +1 ;
    nodes(i)=node;
end

roots = [ xinit, yinit ;
          0.1 , 0.9;
          0.3 , 0.8;
          0.8 , 0.75;
          1.1  , 0.65 ;
          1.1, 0.4 ;
          0.9, 0 ;
          0.7, 0.2 ;
          0.5, 0.4 ;
          0.3, 0.4 ;
          0.2, 0.5 ;
          0.01, 0.6 ;
          0.01, 0.8 ;
          0.7, 0.4;
          0.01, 0.65 ;
          0.01, 0.7 ;
          0.9, 0.8 ;
          1.0, 0.8 ;
          1.1, 0.8 ;
          1, 0.7; 
          1.1, 0.5 ;
          1.1, 0.4 ;
          1.1, 0.3 ;
          0.9, -0.2 
    ];
%plant the roots
nodes(1).x = xinit; 
nodes(1).y = yinit ;
nodes(1).parent = -1; 
nodes(1).direction = +1 ;
tree(1)=nodes(1) ;
r = 1; 
for i=2:num_roots,
    nodes(i).x = roots(i,1);%xinit + r*sin((i-1)*2*pi/(num_roots-1)) ;
    nodes(i).y = roots(i,2);%xinit + r*cos((i-1)*2*pi/(num_roots-1)) ;
    nodes(i).parent = -1; 
    nodes(i).direction = -1 ;
    tree(i)=nodes(i) ;
end


%grow the forest
for i=num_roots+1:num_nodes,
    node = nodes(i);
    indicator = -1 ;
    if  mod(i,100) == 0
        i
    end
   if i>num_nodes/2
      indicator = 0 ;
    end
    nearestNode = findNearestNode(tree, node, indicator);
    node.direction = tree(nearestNode).direction; 
    q_near = [tree(nearestNode).x tree(nearestNode).y] ;
    q_sample = [node.x node.y] ;
    
    %determining a good trajectory
    minimumDistance=9909;
    for j=1:10,
        %If a tree falls in the woods and no one is around to hear it, does it
        %still make a sound? probably not, but it kills schrodinger's cat.
        p = [ ((2*rand)-1)/10;  ((2*rand)-1)/10] ;
        if(node.direction>0)
            [ts, q_new]=ode45(@td2, [0, deltaT], [q_near(1);q_near(2)], options, p);
        else
            [ts, q_new]=ode23s(@td2, [deltaT, 0], [q_near(1);q_near(2)], options, p);
        end
        x= q_new(size(q_new,1),1) ;
        y= q_new(size(q_new,1),2);
        d = norm([(x-q_sample(1)) ; (y-q_sample(2)) ]);
        if(d<minimumDistance),
            minimumDistance=d;
            node.x = x ;
            node.y = y ;
            %node.p = p ;
        end
    end
    node.parent = nearestNode ;
    tree(i) = node; 
   % if node.x > 3 | node.y > 3 | node.x < -3 | node.y < -3 
   %     i = i-1 ;
   % end
end


%this line causes duplicate voronoi vertex to be generated, which later causes some bugs in voronoi figure generation
for i=1:num_nodes,
    if tree(i).x>1.2, 
        tree(i).x=1.2 ;
    end
    if tree(i).x<-0.2, 
        tree(i).x=-0.2;
    end
    if tree(i).y>1.3, tree(i).y=1.3 ;end
    if tree(i).y<-0.3, tree(i).y=-0.3; end
end

hold on
figure(1)
for i=num_roots+1:size(tree,2)
    if tree(i).parent ~= -1,
        if tree(i).direction == +1,
            line([tree(i).x tree(i).x ], [ tree(i).y tree(i).y], 'Marker','.','LineStyle','-', 'Color', 'blue')
             line([tree(i).x tree(tree(i).parent).x ], [ tree(i).y tree(tree(i).parent).y], 'Marker','.','LineStyle','-', 'Color', 'blue' )
        else
            line([tree(i).x tree(i).x ], [ tree(i).y tree(i).y],'Marker','.','LineStyle','-','Color', 'red' )
             line([tree(i).x tree(tree(i).parent).x ], [ tree(i).y tree(tree(i).parent).y],'Marker','.','LineStyle','-','Color', 'red' )
        end
        
    
    %line([tree(i).x tree(i).x+0.01], [ tree(i).y tree(i).y+0.01] )
    end
end


for i=1:num_nodes,
    vx(i) = tree(i).x ;
    vy(i) = tree(i).y ;
end
figure(2)

[vxp, vyp] = voronoi(vx,vy);
plot(vx, vy, 'r+', vxp, vyp, 'b-'); axis equal

figure(3)
vv = [vx;vy]' ;
[v,c] = voronoin(vv);
for i=1:length(c)
    patch(v(c{i},1), v(c{i},2),'w');
end

for i=1:length(c), 
    if tree(i).direction==+1 %all(c{i}~= 1)
        patch(v(c{i},1), v(c{i},2),'y');
    end
end
hold off

%removes unique points
for i=1:num_nodes, 
    for j=1:num_nodes, 
        if i~=j,
            if (tree(i).x==tree(j).x) && (tree(i).y==tree(j).y)
                disp('***')
                i
                j
            end
        end
    end
end

%drawing the tree
%hold on
%for i=1:10,
%   x = rand-0.5 ;
%   y = rand-0.5 ;
%   t = 1 ;
%   [t0,y0] = ode45(@td1,[0, t],[x;y],options);
%   plot(y0(:,1),y0(:,2), 'b');      
   %plot(t0,y0,t1(1:20:end),y1(1:20:end,:),'o')
%end
%for i=1:10,
%   x = rand-0.5 ;
%   y = rand-0.5 ;
%   t = -1e-2;
%   [t1,y1] = ode45(@td1,[t,0],[x;y],options);
%   plot(y1(:,1),y1(:,2), 'r');      
%   %plot(t0,y0,t1(1:20:end),y1(1:20:end,:),'o')
%end
%hold off

%cluster=0; clusterSize=0;
%disp('doing data-mining stuff');
%if config.generateRandomNodes==0,
 %   [cluster clusterSize] = datamining(config);
%end
%disp('Generating node sequence');
%nodes = generateNewNode(config, cluster, clusterSize);
%disp('All-done, now generating RRT')
%tree = generateTree(config, nodes);
%disp('finished, drawing the tree');
%drawTree(tree);


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
