function [tree, flag] = generateTree2( config, tree , xseed, yseed, iteration)
%generateTree The RRT Algorithm
%   The RRT Algorithm

for i=1:iteration
    node.x = xseed + randn/10;
    node.y = yseed + randn/10;
    node.parent = -1 ;
    nearestNode = findNearestNode(tree, node);
    q_near = [tree(nearestNode).x tree(nearestNode).y] ;
    q_sample = [node.x node.y] ;
    %plannar constraints:
    %q_new = q_near + (config.segment/norm(q_near-q_sample)).*(q_sample-q_near) ;
    
    %differential constraints:
    %choosing the best trajectory from q_near towards q_sample
    minimumDistance=config.MAX;
    for j=1:10,
        p = [ randn/10; randn/10 ] ;
        [ts, q_new]=ode45(@josephson, [0,config.deltaT], [q_near(1);q_near(2)], [], p);
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
    tree(size(tree,2)+1) = node; 
end
 
end

