%Josephson performance study
%MO-RRT run until finding 3,00 goal tests
% we report
%   final time,
%   total goal samples
tic
%%
clc;
clear;
config=generateConfig();
%%
goal_mu = [0, 0];
goal_sigma = [0.2 0; 0 0.2];
goal_samples = [mvnrnd(goal_mu, goal_sigma, config.iterations)];

%%
% RRT loop
dummy.x = 0;
dummy.y= 0 ;
dummy.parent=0;
tree(config.iterations)=dummy;
root.x = config.initX ;
root.y= config.initY ;
root.parent=0;
tree(1)=root;
%%

for i=2:config.iterations,
    i
    sample = goal_samples(i,:);
    node.x=sample(1,1); 
    node.y=sample(1,2); 
    node.parent=-1;
        
    nearestNode = findNearestNode(tree, node, i);
    q_near = [tree(nearestNode).x tree(nearestNode).y] ;
    q_sample = [node.x node.y] ;
    minimumDistance=config.MAX;
    for j=1:3,
        p = [ ((2*rand)-1)/10;  ((2*rand)-1)/10] ;
        [ts, q_new]=ode45(@josephson, [0,config.deltaT], [q_near(1);q_near(2)], [], p);
        x= q_new(size(q_new,1),1) ;
        y= q_new(size(q_new,1),2);
        d = norm([(x-q_sample(1)) ; (y-q_sample(2)) ]);
        if(d<minimumDistance),
            minimumDistance=d;
            node.x = x ;
            node.y = y ;
        end
    end
    node.parent = nearestNode ;
    tree(i) = node; 
end
figure(2)
%finalization
disp('finished, drawing the tree');

toc

%scatter(goal_samples(:,1),goal_samples(:,2),10,'.')
drawTree(tree);

%%
total_goal_test=0;
for i=1:config.iterations, 
    x = tree(i).x ;
    y = tree(i).y; 
    d=sqrt(x^2+y^2);
    if(d < 0.5),
        %check if this is a leaf node
        this_is_a_leaf=0;
        for j=1:config.iterations, 
            if(tree(j).parent==i),
                this_is_a_leaf=1;
            end
        end
        if(this_is_a_leaf==0),
            total_goal_test=total_goal_test+1;
        end
    end
end
total_goal_test

%%