%Josephson performance study
%Monte Carlo run until 3,000 iterations
% we report
%   final time,
%   total goal samples
tic
%%
clc;
clear;
config=generateConfig();
%%
dummy.x = 0;
dummy.y= 0 ;
mc(config.iterations)=dummy;
root.x = config.initX ;
root.y= config.initY ;
mc(1)=root;
%%
mcSim=1000;
mcSamplePerSimulation = 250;
iteration=1;
total_goal_samples=0;

for i=1:mcSim,
    iteration
    node.x=config.initX;
    node.y=config.initY;
    for j=1:mcSamplePerSimulation,
        p = [ ((2*rand)-1)/5;  ((2*rand)-1)/5] ;  %input
        [ts, node_new_data]=ode45(@josephson, [0,config.deltaT], [node.x;node.y], [], p);
        node.x =  node_new_data(size(node_new_data,1),1);
        node.y =  node_new_data(size(node_new_data,1),2);
        mc(iteration+1)=node;
        iteration=iteration+1;
    end
    x=node.x;
    y=node.y;
    d=sqrt(x^2+y^2);
    if(d < 0.5),
        total_goal_samples=total_goal_samples+1;
    end
        
end

%finalization
disp('finished, drawing the tree');

toc
total_goal_samples
%scatter(goal_samples(:,1),goal_samples(:,2),10,'.')
%drawTree(tree);

%%
for i=1:mcSim,
    for j=2:mcSamplePerSimulation-1,
        c=(i-1)*mcSamplePerSimulation + j+1;
        line([mc(c).x mc(c-1).x ], [ mc(c).y mc(c-1).y] )
    end
end

%%
