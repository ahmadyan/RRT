%this is the alternative script for RRT journal, jan-2014
%the rrt# algorithm contains a supervised learning step and consists of the
%following sections:
% 1. initial learning phase to identify goal distributions G
% 2. the RRT sample generation from G and H distribution
% 3. RRT growth (standard RRT algorithm)
% 4. extracted reached distribution H from the RRT

clc;
clear;

global variation0
global variation1
config=generateConfig();
sig = config.sigma * ones( 1, config.dim);
var = diag( sig, 0); %=[config.sigma config.csigma; config.csigma config.sigma];

% step 1: identifying goal distribution G
disp('Identifying the goal regions...');

disp('doing data-mining stuff');
load ring_training_sample; %[training_sample] = training(config);
training_sample = training_sample' ;
[label, model, L] = vbgm(training_sample, config.maxcluster); 
[training_dim,training_size] = size(training_sample);
figure(1)
%spread(training_sample, label)

[~,label(:)] = max(model.R,[],2);
goal_index = unique(label);
goal_count = size(goal_index,2);
for i=1:goal_count,
    goal_dist_mu(i,:) = model.m(:,goal_index(i));
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

%initially (for the first 10 sample) we use the same dist for goal and rrt
rrt_dist_mu = goal_dist_mu;
rrt_dist_sigma = goal_dist_sigma;
rrt_dist_p = goal_dist_p;
rrt_dist = gmdistribution(rrt_dist_mu,rrt_dist_sigma,rrt_dist_p);


disp('finished identifying the goal regions.');

%step 1.5: the main simulation path
%We run an initial M.C. random walk without any source of variation in the
%system, this helps us identify the main simulation path of the circuit. 
%We later use this simulation path to determine the signal deviation (i.e.
% |signal_envlope - true_signal|

disp('Initial MC simulation');


%initialing the data, by pre-allocating we can generate data much faster
dummyNode.y = zeros(1, config.dim);
dummyNode.parent = -1; 
dummyNode.t=0;
mc(config.explore_step)=dummyNode;

%initial state
mc_root.y = config.init ;
mc_root.parent=0;
mc_root.t=0;
mc(1)=mc_root;
mc_sample=mc_root.y;

h=config.tmax/config.explore_step;
tmp=mc_root.y; %zeros(1,config.dim); 
t=0;
variation0 = 0; %2*(rand-0.5)/10;
variation1 = 0; %2*(rand-0.5)/20;
   
for i=2:config.explore_step,
    str = sprintf(' Monte Carlo, iter= %d , step=%f.', i, h);
    disp(str);
	tmp = ring(tmp, t, h);
	t=t+h;
    node.y = tmp;
    node.parent = i-1;
    node.t = t;
    mc(i) = node; 
    mc_sample=[ mc_sample; node.y]; 
   
end


%step 2: Growing the RRT in the state space
%initialing the data, by pre-allocating we can generate data much faster
dummyNode.y = zeros(1, config.dim);
dummyNode.parent = -1; 
dummyNode.t=0;
tree(config.iterations)=dummyNode;
% RRT loop
root.y = config.init ;
root.parent=0;
root.t=0;
tree(1)=root;
rrt_sample=root.y;


inputSeq0 = zeros(1, config.iterations);
inputSeq1 = zeros(1, config.iterations);



timeEnvlope=1e-5; 

for i=2:config.iterations,
    config.i = i;
    str = sprintf(' RRT, iter= %d     timeEnvlope=%f', i, timeEnvlope);
    disp(str);
	
    %updating the variance matrices
    sigma = config.final_variance + config.init_variance * exp( -config.variance_cooling_rate * i/config.iterations);
    sig = config.sigma * ones( 1, config.dim );
    var = diag( sig, 0);
    for j=1:goal_count,
        if(j==1)
            goal_dist_sigma = cat( 3, var );
        else
            goal_dist_sigma = cat( 3, goal_dist_sigma, var );
        end
    end
    
    %computing the mixture distribution
    %dist_mu=[goal_dist_mu; rrt_dist_mu];
    %goal_dist_sigma
    %rrt_dist_sigma
    %dist_sigma=cat(3, goal_dist_sigma, rrt_dist_sigma);
    %rrt_dist_size = size(rrt_dist_mu,1);
    %goal_dist_size = size(goal_dist_mu,1);
    %dist_p=[config.alpha*goal_dist_p  (1-config.alpha)*rrt_dist_p];
    %dist=gmdistribution(dist_mu,dist_sigma,dist_p);
    %sample=random(dist, 1); 
    %sample
    
    node.y = zeros(1, config.dim);
    node.parent=-1;
    node.t = config.tmax;
    
    p=rand;
    % generate a new sample point
    if(p<=0.0), %classic RRT algorithm
        for j=1:config.dim,
            node.y(j) = (config.max(j)-config.min(j))*rand + config.min(j);
        end
        node.parent  = -1 ;
    elseif(p<=0.5), %goal-oriented RRT algorithm
        sample = random(goal_dist,1);
        node.y = sample;
        node.parent  = -1 ;
    else
        sample = random(rrt_dist,1);
        node.y = sample;
         for j=5:config.dim,  %the first four dimensions are being determined from the mixture dist
            node.y(j) = (config.max(j)-config.min(j))*rand + config.min(j);
        end
        node.parent  = -1 ;
    end
    
    
    %generate a random sample
    if(i>100),
        for j=5:config.dim,     %the first four dimensions are being determined from the mixture dist.
            node.y(j) = (config.max(j)-config.min(j))*rand + config.min(j);
        end
    else
        for j=1:config.dim,
            node.y(j) = (config.max(j)-config.min(j))*rand + config.min(j);
        end
    end
    
    
    
    %searching for the nearest node in the RRT from q_sample
    
    %setting the time-envlope, if we propagated the envlop until tmax, then
    %discard the timeenvlope and proceed like normal rrt
    tenv=timeEnvlope;
    if(timeEnvlope>= config.tmax),
        tenv=-1;
    end
    nearestNode = findNearestNode(tree, node, config, tenv);
    q_near = tree(nearestNode).y;
    q_sample = node.y ;
    minimumDistance=config.MAX;
    
    %simulating the circuit for q_near
%{  
    %shooting  for choosing an optimum path from q_near towards q_sample
    for j=1:10,
        variation0 = 2*(rand-0.5)/10;
        variation1 = 2*(rand-0.5)/20;
        q_new = ring(q_near, tree(nearestNode).t, config.deltaT);
        x= q_new(size(q_new,1),1) ;
        y= q_new(size(q_new,1),2);
        d = norm([(x-q_sample(1)) ; (y-q_sample(2)) ]);
        if(d<minimumDistance),
           minimumDistance=d;
            node.x = x ;
            node.y = y ;
        end
    end.
%}
    
	variation0 = 2*(rand-0.5)/10;
	variation1 = 2*(rand-0.5)/20;
    inputSeq0(i) = variation0; 
    inputSeq1(i) = variation1;
    
	q_new = ring(q_near, tree(nearestNode).t, config.deltaT);
    
    %Adding the new node the the RRT
    node.y = q_new ;
    node.parent = nearestNode ;
    node.t = tree(nearestNode).t + config.deltaT ;
    tree(i) = node; 
    rrt_sample=[ rrt_sample; node.y]; 
    
    if ( (node.t>timeEnvlope) && (timeEnvlope<=config.tmax)),
        timeEnvlope=node.t;
    end
     
    
    %Computing the distribution of the RRT nodes and the mixture dist.
    if(i>100),
        %variational inference on rrt
        test = rrt_sample(:, [1:4]);
        [label, model, L] = vbgm(test', config.maxcluster); 
        [rrt_sample_dim,rrt_sample_size] = size(rrt_sample');
        %spread(rrt_sample', label)
        [~,label(:)] = max(model.R,[],2);
        rrt_index = unique(label);
        rrt_count = size(rrt_index,2);
        clear rrt_dist_mu; 
        clear rrt_dist_sigma;
        clear rrt_dist_count;
        sig = config.sigma * ones( 1, size(test, 2));
        var = diag( sig, 0);
        for j=1:rrt_count,
            rrt_dist_mu(j, :) = model.m(:,rrt_index(j));
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
drawTree(tree, 2, 1);