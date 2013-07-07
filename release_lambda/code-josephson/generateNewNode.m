function nodes = generateNewNode( config, cluster, clusterSize )
%generateNewNode this function will generate new node for RRT algorithm
% This function will generate new node for RRT algorithm. The new node
% can either be generated radnomly (standard-RRT algorithm) or generated
% from data mining feedback (guided-RRT algorithm).

if config.generateRandomNodes==1,
    for i=1:config.iterations,
        % randomly pick configuration
        node.x       = (config.xmax-config.xmin)*rand + config.xmin;
        node.y       = (config.ymax-config.ymin)*rand + config.ymin;
        node.parent  = -1 ;
        nodes(i)=node;
    end
    
else
    k=1;
    n=100;
    for i=1:config.numberOfClusters,
       fprintf('Analyzing cluster %d\n', i);
        for j=1:n,
            grid = ceil( clusterSize(i)*rand() ) ;
            if grid==0,
                grid=1;
            end
            c = cluster( i, grid ) ;
            x = floor(c/40) - 8;
            y = mod(c,40) - 8;
           
            %todo: scale x and y to system's space
            node.x = x/20 + randn/10 ;
            node.y = y/20 + randn/10 ;
            node.parent = -1 ;
            nodes(k)=node ;
            k=k+1;
        end
        if(i>3)
            n=n*2;
        end
        n = n*2;
    end
    
end

end

