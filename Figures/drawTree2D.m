function drawTree2D( data , ix, iy, xname, yname, xaxismin, xaxismax, yaxismin, yaxismax)
    xmax=-999;
    xmin=999;
    
    ymax=-999;
    ymin=999;
    for i=1:size(data, 2),
        if(data( ix, i) < xmin),
            xmin=data(ix, i);
        end
        
        if(data( ix, i) > xmax),
            xmax=data(ix, i);
        end
        
        if(data( iy, i) < ymin),
            ymin=data(iy, i);
        end
        
        if(data( iy, i) > ymax),
            ymax=data(iy, i);
        end
    end
    
    hold on
    axis([xaxismin xaxismax yaxismin yaxismax])
    set(gca,'FontSize',16,'fontWeight','demi')
    set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

    xlabel('time (s)'); 
    xlabel(xname)
    ylabel(yname)
    
    for i=2:size(data, 2),
        p=data(2, i);
        px=data(ix, p);
        py=data(iy, p);
        x=data(ix,i);
        y=data(iy,i);
        line([px, x], [py, y]);
    end
    hold off
end
