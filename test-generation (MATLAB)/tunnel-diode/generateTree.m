function [tree, flag] = generateTree( config, nodes )
%generateTree The RRT Algorithm
%   The RRT Algorithm

nodes(1).x = config.initX ;
nodes(1).y = config.initY ;
root = nodes(1);
root.parent=0;
root.x = config.initX ;
root.y= config.initY ;
tree(1)=root;

for i=2:size(nodes,2),
    if(mod(i,1000)==0) fprintf('*');
    end
    node = nodes(i) ;
    nearestNode = findNearestNode(tree, node);
    q_near = [tree(nearestNode).x tree(nearestNode).y] ;
    q_sample = [node.x node.y] ;
    %plannar constraints:
    %q_new = q_near + (config.segment/norm(q_near-q_sample)).*(q_sample-q_near) ;
    
    %differential constraints:
    %choosing the best trajectory from q_near towards q_sample
    minimumDistance=config.MAX;
    for j=1:10,
        p = [ ((2*rand)-1)/10;  ((2*rand)-1)/10] ;
        [ts, q_new]=ode45(@tunnel_diode, [0,config.deltaT], [q_near(1);q_near(2)], [], p);
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
    
    %p = [ ((2*rand)-1)/10;  ((2*rand)-1)/10] ;
    %   [ts, q_new]=ode45(@tunnel_diode, [0,config.deltaT], [q_near(1);q_near(2)], [], p);
    %node.x = q_new(size(q_new,1),1) ;
    %node.y = q_new(size(q_new,1),2);
    
    node.parent = nearestNode ;
    tree(i) = node; 
end
 
end

