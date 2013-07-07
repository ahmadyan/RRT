%this is the alternative script for RRT journal, jan-2013
%the rrt# algorithm contains a supervised learning step and consists of the
%following sections:
% 1. initial learning phase to identify goal distributions G
% 2. the RRT sample generation from G and H distribution
% 3. RRT growth (standard RRT algorithm)
% 4. extracted reached distribution H from the RRT
%clc;
%clear;
%config=generateConfig();

var=[config.sigma config.csigma; config.csigma config.sigma];

% step 1: identifying goal distribution G

disp('doing data-mining stuff');
[training_sample] = datamining2(config);
training_sample = training_sample' ;
[label, model, L] = vbgm(training_sample, config.maxcluster); 
[training_dim,training_size] = size(training_sample);
figure(1)
spread(training_sample, label)

[~,label(:)] = max(model.R,[],2);
goal_index = unique(label);
goal_count = size(goal_index,2);
for i=1:goal_count,
    goal_dist_mu(i,1) = model.m(1,goal_index(i));
    goal_dist_mu(i,2) = model.m(2,goal_index(i));
    goal_dist_count(i)=sum(label==goal_index(i));  %number of samples in cluster
    if(i==1)
        goal_dist_sigma = cat( 3, var );
    else
        goal_dist_sigma = cat( 3, goal_dist_sigma, var );
    end
end
%exact answer should (0.063, 0.758) and (0.884, 0.21)

%Guassian Mixture Distribution for Sampling
goal_dist_p =  goal_dist_count/config.learningSamples;
goal_dist = gmdistribution(goal_dist_mu,goal_dist_sigma,goal_dist_p);
%goal_samples = random(goal_dist,config.iterations);

%drawing goal-samples
%ezcontour(@(x,y)pdf(obj,[x y]),[-10 10],[-10 10])
%hold on
%scatter(G_samples(:,1),G_samples(:,2),10,'.')
%ezsurf(@(x,y)pdf(goal_dist,[x y]),[-0.4 1.4],[-0.4 1.4])
% RRT loop
root.x = config.initX ;
root.y= config.initY ;
root.parent=0;
tree(1)=root;

rrt_sample=[root.x root.y]; 

%initially (for the first 10 sample) we use the same dist for goal and rrt
rrt_dist_mu = goal_dist_mu;
rrt_dist_sigma = goal_dist_sigma;
rrt_dist_p = goal_dist_p;
rrt_dist = gmdistribution(rrt_dist_mu,rrt_dist_sigma,rrt_dist_p);

for i=2:config.iterations,
    
    %updating the variance matrices
    sigma = config.final_variance + config.init_variance * exp( -config.variance_cooling_rate * i/config.iterations)
    var=[sigma config.csigma; config.csigma sigma];
    for j=1:goal_count,
    if(j==1)
        goal_dist_sigma = cat( 3, var );
    else
        goal_dist_sigma = cat( 3, goal_dist_sigma, var );
    end
    end
    
    
    fprintf('%d\n', i);
    
    %computing the mixture distribution
    dist_mu=[goal_dist_mu; rrt_dist_mu];
    dist_sigma=cat(3, goal_dist_sigma, rrt_dist_sigma);
    rrt_dist_size = size(rrt_dist_mu,1);
    goal_dist_size = size(goal_dist_mu,1);
    dist_p=[config.alpha*goal_dist_p  (1-config.alpha)*rrt_dist_p];
    dist=gmdistribution(dist_mu,dist_sigma,dist_p);
    
    sample=random(dist, 1); 
    node.x=sample(1,1); 
    node.y=sample(1,2); 
    node.parent=-1;
    %{
    p=rand;
    % generate a new sample point
    if(p<=0.033), %classic RRT algorithm
        node.x       = (config.xmax-config.xmin)*rand + config.xmin;
        node.y       = (config.ymax-config.ymin)*rand + config.ymin;
        node.parent  = -1 ;
    elseif(p<=0.66), %goal-oriented RRT algorithm
        sample = random(goal_dist,1);
        node.x = sample(1, 1);
        node.y = sample(1, 2);
        node.parent  = -1 ;
    else
        sample = random(rrt_dist,1);
        node.x = sample(1, 1);
        node.y = sample(1, 2);
        node.parent  = -1 ;
    end
    %}  
        
    
    
    nearestNode = findNearestNode(tree, node);
    q_near = [tree(nearestNode).x tree(nearestNode).y] ;
    q_sample = [node.x node.y] ;
    minimumDistance=config.MAX;
    for j=1:10,
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
    rrt_sample=[ rrt_sample; node.x node.y]; 
     
    if(i>10),
        %variational inference on rrt
        [label, model, L] = vbgm(rrt_sample', config.maxcluster); 
        [rrt_sample_dim,rrt_sample_size] = size(rrt_sample');
        %spread(rrt_sample', label)
        [~,label(:)] = max(model.R,[],2);
        rrt_index = unique(label);
        rrt_count = size(rrt_index,2);
        clear rrt_dist_mu; 
        clear rrt_dist_sigma;
        clear rrt_dist_count;
        for j=1:rrt_count,
            rrt_dist_mu(j, 1) = model.m(1,rrt_index(j));
            rrt_dist_mu(j, 2) = model.m(2,rrt_index(j));
            rrt_dist_count(j)=sum(label==rrt_index(j));  %number of samples in cluster
            if(j==1)
                rrt_dist_sigma = cat( 3, var );
            else
                rrt_dist_sigma = cat( 3, rrt_dist_sigma, var );
            end
        end
        rrt_dist_p =  rrt_dist_count /i;
        rrt_dist = gmdistribution(rrt_dist_mu,rrt_dist_sigma,rrt_dist_p);
    end    
end
figure(2)
%finalization
disp('finished, drawing the tree');
drawTree(tree);