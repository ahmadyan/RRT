function core = parser( config , core )
%PARSER Summary of this function goes here
%   Detailed explanation goes here

%% Initialize
numElem=0;  %Number of passive elements.
numV=0;     %Number of independent voltage sources
numO=0;     %Number of op amps
numI=0;     %Number of independent current sources
numD=0;

numE=0;
numF=0;
numG=0;
numH=0;
numQ=0;

numNode=0;  %Number of nodes, not including ground (node 0).
node1 =0;
node2 =0;

components = java.util.Hashtable ;
nodes = java.util.Hashtable ;
components.clear;
nodes.clear;
nodes.put ('0', 0); %Adding ground
numNode = 1 ;

Element=[];
Vsource=[];
Isource=[];
Esource=[];
Fsource=[];
Gsource=[];
Hsource=[];
Diodes=[];
tran=[];
vNodeset=[];
Nodeset=[];

[Name arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9]= textread(config.fileName,'%s %s %s %s %s %s %s %s %s %s') ; 


%% Parsing spice netlist
for i=1:length(Name),
    if(config.verbose==1) disp(Name{i}) ; end
    switch(Name{i}(1)),
        case {'R','L','C','r','l','c'},
            components.put(Name{i} , numElem);
            numElem=numElem+1;
            Element(numElem).Name=Name{i};
            Element(numElem).Node1=getNodeNumber( nodes, arg1{i} );
            Element(numElem).Node2=getNodeNumber( nodes, arg2{i} );
            
            try
                Element(numElem).Value=str2double(arg3{i});
            catch
                Element(numElem).Value=nan;
            end
        case {'V', 'v'},
            numV=numV+1;
            Vsource(numV).Name=Name{i};
            Vsource(numV).Node1=getNodeNumber( nodes, arg1{i} );
            Vsource(numV).Node2=getNodeNumber( nodes, arg2{i} );
            
            if(arg3{i}(1)=='a')
                Vsource(numV).type = 1 ; %VAC
                Vsource(numV).peak = str2double (substr('ac=0.0001', 4,10) ) ;
                disp(sprintf('[warning] ac-source detected %s ', Name{i},  arg3{i})); 
            elseif(arg3{i}(1)=='S')
                Vsource(numV).type = 2 ; %VSIN
                Vsource(numV).offset=str2double(arg4{i});
                Vsource(numV).peak=str2double(arg5{i});
                Vsource(numV).freq=str2double(arg6{i});
            else
                Vsource(numV).type = 0 ; %VDC
                
            try
                Vsource(numV).Value=str2double(arg3{i});
            catch
                Vsource(numV).Value=nan;
            end
            
            end
            
        %case {'O', 'o' },
        %    numO=numO+1;
        %    Opamp(numO).Name=Name{i};
        %    Opamp(numO).Node1=str2double(arg1{i});
        %    Opamp(numO).Node2=str2double(arg2{i});
        %    Opamp(numO).Node3=str2double(arg3{i});
        
        case {'I', 'i'},
            numI=numI+1;
            Isource(numI).Name=Name{i};
            Isource(numI).Node1=getNodeNumber( nodes, arg1{i} );
            Isource(numI).Node2=getNodeNumber( nodes, arg2{i} );
            try
                Isource(numI).Value=str2double(arg3{i});
            catch
                Isource(numI).Value=nan;
            end
            
        %Linear Voltage-Controlled Voltage Sources:	EXXXXXXX N+ N- NC+ NC- VALUE
        %example: e1 c 0 a b 4
        case {'E', 'e'},
            numE=numE+1;
            Esource(numE).Name = Name{i};
            Esource(numE).Node1 = getNodeNumber( nodes, arg1{i} );
            Esource(numE).Node2 = getNodeNumber( nodes, arg2{i} );
            Esource(numE).Node3 = getNodeNumber( nodes, arg3{i} );
            Esource(numE).Node4 = getNodeNumber( nodes, arg4{i} );
            try
            	Esource(numE).Value=str2double(arg5{i});
            catch
                Esource(numE).Value=nan;
            end
            
        
        %Linear Current-Controlled Current Sources:	FXXXXXXX N+ N- VNAM VALUE
        %example: f1 0 e v1 8
        case {'F', 'f'},
            numF=numF+1;
            Fsource(numF).Name = Name{i};
            Fsource(numF).Node1 = getNodeNumber( nodes, arg1{i} );
            Fsource(numF).Node2 = getNodeNumber( nodes, arg2{i} );
            Fsource(numF).Dependant = arg3{i};
            try
            	Fsource(numF).Value=str2double(arg4{i});
            catch
                Fsource(numF).Value=nan;
            end
            
        %Linear Voltage-Controlled Current Sources:	GXXXXXXX N+ N- NC+ NC- VALUE
        %example: g9 g h e f 13
        case {'G', 'g'},
            numG=numG+1;
            Gsource(numG).Name = Name{i};
            Gsource(numG).Node1 = getNodeNumber( nodes, arg1{i} );
            Gsource(numG).Node2 = getNodeNumber( nodes, arg2{i} );
            Gsource(numG).Node3 = getNodeNumber( nodes, arg3{i} );
            Gsource(numG).Node4 = getNodeNumber( nodes, arg4{i} );
            try
            	Gsource(numG).Value=str2double(arg5{i});
            catch
                Gsource(numG).Value=nan;
            end
            
          
        %Linear Current-Controlled Voltage Sources:	HXXXXXXX N+ N- VNAM VALUE
        %example: h1 g 0 v1 11
        case {'H', 'h'},
            numH=numH+1;
            Hsource(numH).Name = Name{i};
            Hsource(numH).Node1 = getNodeNumber( nodes, arg1{i} );
            Hsource(numH).Node2 = getNodeNumber( nodes, arg2{i} );
            Hsource(numH).Dependant = arg3{i};
            try
            	Hsource(numH).Value=str2double(arg4{i});
            catch
                Hsource(numH).Value=nan;
            end
      
            
        case {'D', 'd'},
            numD=numD+1;
            Diodes(numD).Name = Name{i};
            Diodes(numD).Node1 = getNodeNumber( nodes, arg1{i} );
            Diodes(numD).Node2 = getNodeNumber( nodes, arg2{i} );
            
            
            components.put(Name{i} , numElem);
            numElem=numElem+1;
            Element(numElem).Name=Name{i};
            Element(numElem).Node1=getNodeNumber( nodes, arg1{i} );
            Element(numElem).Node2=getNodeNumber( nodes, arg2{i} );
            
            try
                Element(numElem).Value=str2double(arg3{i});
            catch
                Element(numElem).Value=nan;
            end
            
            
        case {'Q', 'q'},
            numQ=numQ+1;
            tran(numQ).Name = Name{i};
            tran(numQ).Node1 = getNodeNumber( nodes, arg1{i} );
            tran(numQ).Node2 = getNodeNumber( nodes, arg2{i} );
            tran(numQ).Node3 = getNodeNumber( nodes, arg3{i} );
            tran(numQ).Type = arg4{i};
            
            
            components.put(Name{i} , numElem);
            numElem=numElem+1;
            Element(numElem).Name=Name{i};
            Element(numElem).Node1=getNodeNumber( nodes, arg1{i} );
            Element(numElem).Node2=getNodeNumber( nodes, arg2{i} );
            Element(numElem).Node3=getNodeNumber( nodes, arg3{i} );
            Element(numElem).Type= arg4{i} ;
            
    end
    %numNode=max(str2num(arg1{i}),max(str2num(arg2{i}),numNode));
end
% Parsing completed

%% Parsing Spice commands
vNodeset = zeros(size(nodes)-1,1);
iNodeset = zeros(numElem+numV+numI+numE+numF+numG+numH+numD+numQ,1);

for i=1:length(Name),
    switch(Name{i}(1)),
        case {'.'},
            if(strncmpi(Name{i}, '.TRAN',10)==1)
                core.dt=0.1e-6;
                core.t0=0;
                core.tf=200e-6;
            elseif(strncmpi(Name{i}, '.nodeset',10)==1)
                %arg1
                r = strread(regexprep(arg1{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                %arg2
                r = strread(regexprep(arg2{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                
                
                r = strread(regexprep(arg3{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                
                r = strread(regexprep(arg4{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                
                r = strread(regexprep(arg5{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                 r = strread(regexprep(arg6{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                 r = strread(regexprep(arg7{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                 r = strread(regexprep(arg8{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                 r = strread(regexprep(arg9{i},'[=() ]',char(1)),'%s','delimiter',char(1));
                if(length(r)~=0)
                if (strncmpi(r(1), 'v',2)==1)
                    vNodeset( getNodeNumber(nodes, char( r(2)) ) ) = str2double ( char( r(4) ) ) ;
                else
                    for j=1:length(Name),
                        if (strncmpi(Name{j}, r(2),10)==1)
                            iNodeset(j) = str2double ( char( r(4) ) ) ;
                        end
                    end
                end
                end
                
                
            end
    end
end



%% Post processing 
    %processing current-controlled voltage source 
    dependantId = 0;
    for i=1:numH,
        DependentSourceName = Hsource(i).Dependant ;
        j=0;
        for j=1:numV,
            if(strncmpi(Vsource(i).Name, Hsource(i).Dependant, 10)==1) dependantId = j ;
            end
        end
        if(dependantId==0)
            disp(sprintf('[error] netlist error. the dependant source %s for source %s not found', Hsource(i).Dependant, Hsource(i).Name ));
        end
        Hsource(i).Node3 = Vsource(dependantId).Node1 ;
        Hsource(i).Node4 = Vsource(dependantId).Node2 ;
        Hsource(i).DependantID = dependantId ;
    end
    
    
    
    %processing current-controlled current source 
    dependantId = 0;
    for i=1:numF,
        DependentSourceName = Fsource(i).Dependant ;
        j=0;
        for j=1:numV,
            if(strncmpi(Vsource(i).Name, Fsource(i).Dependant, 10)==1) dependantId = j ;
            end
        end
        if(dependantId==0)
            disp(sprintf('[error] netlist error. the dependant source %s for source %s not found', Hsource(i).Dependant, Hsource(i).Name ));
        end
        Fsource(i).Node3 = Vsource(dependantId).Node1 ;
        Fsource(i).Node4 = Vsource(dependantId).Node2 ;
        Fsource(i).DependantID = dependantId ;
    end
% post processing completed.
numNode = size(nodes)-1; %I alaways count one extra!


core.numElem=numElem;  %Number of passive elements.
core.numV=numV;     %Number of independent voltage sources
core.numO=numO;     %Number of op amps
core.numI=numI;     %Number of independent current sources
core.numD=numD;

core.numE=numE;
core.numF=numF;
core.numG=numG;
core.numH=numH;
core.numQ=numQ;

core.numNode=numNode;  %Number of nodes, not including ground (node 0).

core.components = components ;
core.nodes = nodes ;

core.Element=Element;
core.Vsource=Vsource;
core.Isource=Isource;
core.Esource=Esource;
core.Fsource=Fsource;
core.Gsource=Gsource;
core.Hsource=Hsource;
core.Diodes=Diodes;
core.tran=tran;
core.vNodeset=vNodeset;
core.iNodeset=iNodeset;



end

