function drawTrace(data, index, t, yname, tmax, yaxismin, yaxismax)   
    timeEnvlope=0;
    max=-999;
    min=999;
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
    
    
    for i=2:size(data, 2),
        p=data(2, i);
        
        px=data(index, p);
        pt=data(t, p);
        cx=data(index,i);
        ct=data(t,i);
       
        if(tmax>0),
            if(ct<=tmax),
                 line([pt, ct], [px, cx]);
            end
        else
             line([pt, ct], [px, cx]);
            
        end
       
    end
    hold off
end
