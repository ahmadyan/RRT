function drawEye(data, index, t, offset, period, window, yname, yaxismin, yaxismax)   
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
    
    if(offset==-1),
        d0 = data(index, 1);
        i=2;
        if(d0>0.6),
            while data(index, i)>0.5
                i = i+1;
            end
            offset=data(t, i);
        else
            while data(index, i)<0.5
                i = i+1;
            end
            offset=data(t, i);
        end
    end
    offset
    hold on
    %axis([0 2*period yaxismin yaxismax])
    axis([0 window yaxismin yaxismax])
    
    set(gca,'FontSize',16,'fontWeight','demi')
    set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

    xlabel('time (s)'); 
    ylabel(yname)
    
    
    for i=2:size(data, 2),
        p=data(2, i);
            if(p~=1),   %ignore the root
                px=data(index, p);
                pt=data(t, p);
                cx=data(index,i);
                ct=data(t,i);
       
                pt=pt-offset;
                ct=ct-offset;
                clk=ceil(pt/period)-1;
                t0=pt-clk*period;
                t1=ct-clk*period;
        
                %t2=t0+period/2;
                %t3=t1+period/2;
                t2=t0+period;
                t3=t1+period;
        
                
                line([t0, t1], [px, cx]);
                if(t3>t2),
                    line([t2, t3], [px, cx]);
                end
            end
    end
    
     for i=2:size(data, 2),
        p=data(2, i);
            if(p~=1),   %ignore the root
                px=data(index, p);
                pt=data(t, p);
                cx=data(index,i);
                ct=data(t,i);
       
                pt=pt-offset;
                ct=ct-offset;
                clk=ceil(pt/period)-1;
                t0=pt-clk*period;
                t1=ct-clk*period;
        
                %t2=t0+period/2;
                %t3=t1+period/2;
                t2=t0;%+period;
                t3=t1;%+period;
                
                if(t2>period),
                    t2=t2-period;
                end
                if(t3>period),
                    t3=t3-period;
                end
                
                line([t0, t1], [px, cx]);
                if(t3>t2),
                    line([t2, t3], [px, cx]);
                end
            end
     end
    
    hold off
end
