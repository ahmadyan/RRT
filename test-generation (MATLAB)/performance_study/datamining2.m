%alternative datamining for rrt-journal, jan-2013
function [ sample ] = datamining( config )
    sample = zeros(config.learningSamples,2);
    count=1;
    hold on
    %{
    %% previous training code, 
    %% bug: if samples goes outside of the state space "box"
    for i=1:config.learningSamples,
        x = (config.xmax-config.xmin)*rand + config.xmin;
        y = (config.ymax-config.ymin)*rand + config.ymin ;
        t = config.sim_time+3*rand ;
        [ts,ys] = ode45(@josephson,[0,t],[x;y], [], [0;0]); 
        plot(ys(:,1),ys(:,2)); 
    
        %check if it is within the universe:
        sample(count, 1)=ys(size(ys,1),1);
        sample(count, 2)=ys(size(ys,1),2);
        count=count+1;
    end
    %}
    
    while( count <= config.learningSamples ),
        x = (config.xmax-config.xmin)*rand + config.xmin;
        y = (config.ymax-config.ymin)*rand + config.ymin ;
        t = config.sim_time+3*rand ;
        [ts,ys] = ode45(@josephson,[0,t],[x;y], [], [0;0]); 
        plot(ys(:,1),ys(:,2)); 
        
        if ((config.xmin <= ys(size(ys,1),1)) && (ys(size(ys,1),1)<= config.xmax) && (config.ymin <= ys(size(ys,1),2)) && (ys(size(ys,1),2) <= config.ymax)),
            sample(count, 1)=ys(size(ys,1),1);
            sample(count, 2)=ys(size(ys,1),2);
            count=count+1;
        end
    end
    
%{
for i=config.resolution*config.xmin:config.resolution*config.xmax,
    for j=config.resolution*config.ymin:config.resolution*config.ymax,
        [ts,ys] = ode45(config.system,[0,config.sim_time+3*rand],[i*(1/config.resolution);j*(1/config.resolution)]); 
        plot(ys(:,1),ys(:,2)); 
        sample(count, 1)=ys(size(ys,1),1);
        sample(count, 2)=ys(size(ys,1),2);
        count=count+1;
    end
end
vectfield(config.system,config.xmin:.05:config.xmax,config.ymin:.05:config.ymax); 
%}

figure(2)
scatter(sample(:,1), sample(:,2))
hold off
%{

%
%   Grid-based clustering algorithm at DAC-12
%   This has been replaced with Variational Inference algorithm for TCAD
Z=zeros(40,40);
for i=1:size(sample,1),
    x=sample(i,1);
    y=sample(i,2); 
    %Z(4+floor(x*10), 4+floor(y*10)) = Z(4+floor(x*10), 4+floor(y*10))+1;
    Z(8+floor(x*20), 8+floor(y*20)) = Z(8+floor(x*20), 8+floor(y*20))+1;
end

for i=1:40,
    for j=1:40,
        Y(i,j)=Z(i,j);
    end
end


%upward-smoothing
for i=1+1:40-1,
    for j=1+1:40-1,
        x = Y(i,j);
        if(Z(i-1,j-1)<x), 
            Z(i-1,j-1)= ( x + Z(i-1,j-1) )/2;
        end
        
        if(Z(i,j-1)<x), 
            Z(i,j-1)= ( x + Z(i,j-1) )/2;
        end
        
        if(Z(i+1,j-1)<x), 
            Z(i+1,j-1)= ( x + Z(i+1,j-1) )/2;
        end
        
        if(Z(i-1,j)<x), 
            Z(i-1,j)= ( x + Z(i-1,j) )/2;
        end
        
        if(Z(i+1,j)<x), 
            Z(i+1,j)= ( x + Z(i+1,j) )/2;
        end
        
        if(Z(i-1,j+1)<x), 
            Z(i-1,j+1)= ( x + Z(i-1,j+1) )/2;
        end
        
        if(Z(i,j+1)<x), 
            Z(i,j+1)= ( x + Z(i,j+1) )/2;
        end
        
        if(Z(i+1,j+1)<x), 
            Z(i+1,j+1)= ( x + Z(i+1,j+1) )/2;
        end
    end
end

min=1; max=1;

for i=1:40
    for j=1:40
        if Z(i,j)<min, 
            min=Z(i,j);
        end
        if Z(i,j)>max, 
            max=Z(i,j);
        end        
    end
end

%for i=1:20
%    for j=1:20
%        if Z(i,j)<min, 
%            min=Z(i,j);
%        end
%        if Z(i,j)>max, 
%            max=Z(i,j);
%        end
%        
%    end
%end

%grid-based clustering
p = (max-min)/config.numberOfClusters;
clusterSize = zeros(config.numberOfClusters,1);
cluster = zeros(config.numberOfClusters, 40*40);
for i=1:40
    for j=1:40
        c = ceil(Z(i,j)/p);
        if (c==0) 
            c=1;
        end
        cluster ( c, clusterSize(c)+1 ) = 40*i+j;
        clusterSize(c)=clusterSize(c)+1;
    end
end

for i=1:40,
    for j=1:40,
        Y(i,j)=(-Z(j,i))+1;
    end
end



figure(3)
surf(Y);
%}
end

