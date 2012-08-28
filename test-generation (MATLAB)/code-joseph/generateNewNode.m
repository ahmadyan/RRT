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
    n=10;
    for i=1:config.numberOfClusters,
       fprintf('Analyzing cluster %d\n', i);
        for j=1:n,
            grid = ceil( clusterSize(i)*rand() ) ;
            if grid==0,
                grid=1;
            end
            c = cluster( i, grid ) ;
            x = floor(c/10) - 10;
            y = mod(c,10) - 10;
           
            %todo: scale x and y to system's space
            node.x = x + 3*randn ;
            node.y = y + 3*randn ;
            node.parent = -1 ;
            nodes(k)=node ;
            k=k+1;
        end
        if( i>3 ), n = n*2,
        end
        n = n*2;
    end
    
end

end

