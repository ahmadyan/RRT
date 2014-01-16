%alternative datamining for rrt-journal, jan-2013
%updated for a d-dimensional case
function [ sample ] = training( config )
    global variation0
    global variation1
    sample = zeros(config.learningSamples, config.dim);
    count=1;
    variation0=0;
	variation1=0;
    
    while( count <= config.learningSamples ),
        y = zeros(1, config.dim);
        for j=1:config.dim,
            y(j) = (config.max(j)-config.min(j))*rand + config.min(j);
        end
        count
        t = config.tmax*rand + config.sim_time ;
        ys = ring(y, 0, t);
        
        %plot(ys(:,1),ys(:,2)); 
        
        invalidsample=false;
        for j=1:config.dim,
            if(y(j)>config.max(j) || y(j)<config.min(j) ),
                invalidsample=true;
            end
        end
        if(~invalidsample),
            sample(count, :)=ys;
            count=count+1;
        end
    end
    
scatter(sample(:,1), sample(:,2))


end

