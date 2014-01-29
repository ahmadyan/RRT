function drawTrace(data, index, t, yname)
    
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
    
    
    figure(1)
    hold on
    axis([0 timeEnvlope min max])
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
        line([pt, ct], [px, cx]);
       
    end
    hold off
end
