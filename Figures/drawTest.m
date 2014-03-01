%   Draw a test sequence from sample leaf to the root of the tree
%
%
function drawTest(data, index, t, leaf, yname, tmax, yaxismin, yaxismax, t0)   
    timeEnvlope=0;
    max=-999;
    min=999;
    
    if(t0>0),
    tt=0;
    i=1;
    while(tt<t0),
        if( data( t, i) > tt),
            leaf=i;
            tt=data(t,i);
        end
        i=i+1;
    end
    
    end
    
    for i=1:size(data, 2),
        if(data( t, i) > timeEnvlope),
            timeEnvlope=data(t, i);
        end
        
        if(data( index, i) < min),
            min=data(index, i);
        end
        
        if(data( index, i) > max),
            max=data(index, i);
        end
        
        
    end
    
    
    hold on
    if(tmax>0),
        axis([0 tmax yaxismin yaxismax])
    else
        axis([0 timeEnvlope yaxismin yaxismax])
    end
    
    set(gca,'FontSize',16,'fontWeight','demi')
    set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

    xlabel('time (s)'); 
    ylabel(yname)
    
    current=leaf;
    parent=data(2, leaf);
    while parent~=0,
        currentX=data(index, current);
        currentT=data(t,current);
        parentX=data(index,parent);
        parentT=data(t,parent);
        
        
        if(tmax>0),
            if(ct<=tmax),
                 line([parentT, currentT], [parentX, currentX],'color',  'b');
            end
        else
             line([parentT, currentT], [parentX, currentX], 'color', 'b');
            
        end
        
        current=parent;
        parent=data(2,current);
        
    end
    
    hold off
end