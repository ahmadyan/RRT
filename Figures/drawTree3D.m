function drawTree3D( data , ix, iy, it, xname, yname, xaxismin, xaxismax, yaxismin, yaxismax, tmax)
    xmax=-999;
    xmin=999;
    
    ymax=-999;
    ymin=999;
    
    tEnvlope=0;
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
        
        if(data( it, i) > tEnvlope),
            tEnvlope=data(it, i);
        end
    end
    
    if(xaxismin>0),
        x1=xaxismin;
    else
        x1=xmin;
    end
    
    if(xaxismax>0),
        x2=xaxismax;
    else
        x2=xmax;
    end
    
    
    if(yaxismin>0),
        y1=yaxismin;
    else
        y1=ymin;
    end
    
    
    if(yaxismax>0),
        y2=yaxismax;
    else
        y2=ymax;
    end
    
    t1=0;
    if(tmax>0),
        t2=tmax;
    else
        t2=tEnvlope;
    end
    
    
    
    hold on
    
    axis([x1 x2 y1 y2 t1 t2])
    
    set(gca,'FontSize',16,'fontWeight','demi')
    set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
 
    xlabel(xname)
    ylabel(yname)
    zlabel('time (s)');
    
    for i=2:size(data, 2),
        p=data(2, i);
        px=data(ix, p);
        py=data(iy, p);
        pt=data(it, p);
        x=data(ix,i);
        y=data(iy,i);
        t=data(it,i);
        if(t<=t2),
            line([px, x], [py, y], [pt t]);
        end
    end
    hold off
end
