function reportGenerator(config, core, freq, probe)
    
    if(core.mode==0)
    % Fill the V matrix, Holding variable names
        for i=1:core.numNode,
            V{i}=['V(' num2str(i) ')'];
        end

        for i=1:core.numV,
            V{core.numNode+i}=[sprintf('i(%d,%d)',  core.Vsource(i).Node1, core.Vsource(i).Node2)];
        end

        for i=1:core.numE,
            V{i+core.numNode+core.numV}=[sprintf('i(%d,%d)',  core.Esource(i).Node1, core.Esource(i).Node2)] ;
        end

        %Linear Current-Controlled Current Sources
        for i=1:core.numF,
            V{i+core.numNode+core.numE+core.numV} = [sprintf('i(%d,%d)',  core.Fsource(i).Node3, core.Fsource(i).Node4)];
        end

        %Linear Current-Controlled Voltage Sources:
        for i=1:core.numH,
            V{i+core.numNode+core.numE+core.numF+core.numV} = [sprintf('i(%d,%d)',  core.Hsource(i).Node1, core.Hsource(i).Node2)];
        end

        disp(V);

    end
    
    if(core.mode==1)
        
    figure(1);
        %generating a magnitude plot for ac-anlaysis
        %freq&probe comes from ac module
        plot( freq, abs(probe) );
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        figure(2);
        %phase-plot for ac analysis
        plot( freq, angle(probe)*180/pi)
        set(gca,'XScale','log')
    end
    
end