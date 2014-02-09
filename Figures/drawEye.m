function drawEye(data, index, t, period, yname, yaxismin, yaxismax)   
    max=-999;
    min=999;
    for i=1:size(data, 2),
        if(data( index, i) < min),
            min=data(index, i);
        end
        
        if(data( index, i) > max),
            max=data(index, i);
        end
    end
    
    
    hold on
    axis([0 period yaxismin yaxismax])
    
    set(gca,'FontSize',16,'fontWeight','demi')
    set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

    xlabel('time (s)'); 
    ylabel(yname)
    
    
    for i=2:size(data, 2),
        
        p=data(2, i);
        if(p~=1),
        
        px=data(index, p);
        pt=data(t, p);
        cx=data(index,i);
        ct=data(t,i);
       
        
        clk=ceil(pt/period)-1;
        t0=pt-clk*period;
        t1=ct-clk*period;
        
        t2=t0+period/2;
        t3=t1+period/2;
        
        if(t2>period),
            t2 = t2-period;
        end
       
        if(t3>period),
            t3 = t3-period;
        end 
        line([t0, t1], [px, cx]);
        if(t3>t2),
            line([t2, t3], [px, cx]);
        end
       end
    end
    hold off
end
