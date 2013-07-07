%% ECE-552 Spice project / Phase I 
%   Author: Seyed Nematollah Ahmadyan
%   Fall-2011

%   this is a kernel script
tic
clear
disp(sprintf('ECE552-Spice Started')); 
fileName = 'bjt.txt' ;
verbose=0;
symbolic=0;
format long ;
hspiceEnabled=0;
disp(sprintf('[info] parsing input netlist: %s', fileName)); 
if hspiceEnabled, 
    system(sprintf('hspice %s', fileName));
end


bF=100;% forward beta
bR=4; % reverse beta
Va=10;	    % Early voltage
tau = 2 * 10^(-11) ;
Is = 10^(-15);  % saturation current
Vt = 0.02585126075 ;    % 1/40  % thermal voltage
Cj = 10 ^ (-14) ;
Vj = 0.8 ;
fc = 0.5 ;
mj = 0.5 ;

format long

[Name arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9]= textread(fileName,'%s %s %s %s %s %s %s %s %s %s') ; 


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

% initialization completed

%% Parsing spice netlist
for i=1:length(Name),
    if(verbose==1) disp(Name{i}) ; end
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
            if(arg3{i}(1)=='a')
                disp(sprintf('[warning] ignoring %s because of %s', Name{i},  arg3{i})); 
            else
                
            numV=numV+1;
            Vsource(numV).Name=Name{i};
     
            Vsource(numV).Node1=getNodeNumber( nodes, arg1{i} );
            Vsource(numV).Node2=getNodeNumber( nodes, arg2{i} );
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
            if(strncmpi(Name{i}, '.nodeset',10)==1)
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


vNodeset
diff = 1 ;
iteration = 0 ;
%% Main loop for non-linear DC Analysis using Newton-Raphson Method
%while( diff > 10^(-10) )
while (iteration<10)


%% Preallocate all of the cell arrays #################################

%Circuit Matrix (Symbolic and numErical)
if (symbolic==1)
    YS=cell(numNode + numV + numE + numF + numH ,numNode + numV + numE + numF + numH);
end
Y=zeros(numNode + numV + numE + numF + numH ,numNode + numV +  numE + numF + numH);
%Y=sparse(numNode + numV + numE + numF + numH ,numNode + numV + numE + numF + numH);

%Variable Matrix
V=cell(numNode+ numV + numE + numF + numH,1);

if (symbolic==1)
    JS=cell(numNode+ numV + numE + numF + numH,1);
end
J=zeros(numNode+ numV + numE + numF + numH,1);
%J=sparse(numNode+ numV + numE + numF + numH,1);

% Done preallocating cell arrays -------------------------------------

    
%% Fill the G matrix ##################################################
%Initially, make the G Matrix all zeros.
if (symbolic==1)
    [YS{:}]=deal('0');
end

%Now fill the G matrix with conductances from netlist
for i=1:numElem,
    n1=Element(i).Node1;
    n2=Element(i).Node2;
    %Make up a string with the conductance of current element.
    switch(Element(i).Name(1)),
        case {'R' , 'r'},
            g = ['1/' Element(i).Name];
            v = 1/Element(i).Value;
        case {'L', 'l'},
            g = ['1/s/' Element(i).Name];
        case {'C', 'c'},
            g = ['s*' Element(i).Name];
        case {'D', 'd'},
            if(Element(i).Node1 ~= 0)
                v1 = vNodeset( Element(i).Node1  );
            else
                v1 = 0;
            end
            if(Element(i).Node2 ~= 0)
                v2 = vNodeset( Element(i).Node2  );
            else
                v2 = 0;
            end
            
            %diode equations:
            vd=v1-v2;
            id = Is * ( exp(vd/Vt)-1 ); %standard diode equation
            geq = id/Vt ;       % g_eq = dId/dVd = IS/Vt * exp(Vd/vt) = i(v_n)/V_t         
            v = geq;
        case {'Q', 'q'},
            v=0;
            
    end
    
    
    
    %If neither side of the element is connected to ground
    %then subtract it from appropriate location in matrix.
    if (n1~=0) && (n2~=0),
        if (symbolic==1)
            YS{n1,n2}=[ YS{n1,n2} '-' g];
            YS{n2,n1}=[ YS{n2,n1} '-' g];
        end
        Y(n1,n2)= Y(n1,n2) - v;
        Y(n2,n1)= Y(n2,n1) - v;
    end
    
    %If node 1 is connected to ground, add element to diagonal
    %of matrix.
    if (n1~=0),
        if (symbolic==1)
            YS{n1,n1}=[ YS{n1,n1} '+' g];
        end
        Y(n1,n1)= Y(n1,n1) + v;
    end
    %Ditto for node 2.
    if (n2~=0),
        if (symbolic==1)
            YS{n2,n2}=[ YS{n2,n2} '+' g];
        end
        Y(n2,n2)= Y(n2,n2) + v;
    end
end


%The G matrix is finished -------------------------------------------


%% Stamping BJTs

for i=1:numQ,
    c=tran(i).Node1 ;
       b=tran(i).Node2 ;
       e=tran(i).Node3 ;
       if (strncmpi(tran(i).Type, 'npn',4)==1) polarity=1 ;
       else polarity=-1; end
       if(c ~= 0) vc = vNodeset( c  );
       else vc = 0; end
       if(b ~= 0) vb = vNodeset( b  );
       else vb = 0; end
       if(e ~= 0) ve = vNodeset( e  );
       else ve = 0; end
       vbe=polarity*(vb-ve);
       vce=(vc-ve);
       vbc=polarity*(vb-vc);
       
       % gummel-poon transistor model
       % equations derived from http://dev.hypertriton.com/edacious/trunk/generic/npn.c
       Ibf = (Is/bF)*(exp(vbe/Vt)-1);   %Eq. 13
       Ibr = (Is/bR)*(exp(vbc/Vt)-1);   %Eq. 14
       Icc = bF*Ibf-bR*Ibr; %*(1+vce/va)
       gPiF = Ibf/Vt ;                  %Eq. 15, g_ce
       gPiR = Ibr/Vt ;                  %Eq. 16, g_ec
       gmF  = bF*gPiF ;                  % g_ee
       gmR  = bR*gPiR ;                  % g_cc
       go  = Icc/Va;
       Ibf_eq = Ibf - gPiF*vbe;
       Ibr_eq = Ibr - gPiR*vbc;
       Icc_eq = Icc - gmF*vbe + gmR*vbc; %-go*vce
       
       gee=gmF;
       gec=-gPiR;
       gce=-gPiF;
       gcc=gmR;
       ie= polarity*(gee*vbe + gec*vbc) ;
       ic= polarity*(gce*vbe + gcc*vbc) ;
       ib= polarity*((ie+ic));
       
	%disp('----------------BEFORE-----------------')		
      % Y
      % J
     %  disp('----------------AFTER-----------------')
       %Y = StampConductance(    Y,  gPiF,   b,  e   );
       %Y = StampConductance(    Y,  gPiR,   b,  c   );
       %Y = StampConductance(    Y,  go,     e,  c   );
       %Y = StampVCCS(           Y,  gmF,    b,  e,  c,  e 	);
       %Y = StampVCCS(           Y,  gmR,    b,  c,  e,  c   );
       %J = StampCurrentSource(  J,  Ibf_eq, e,  b   );
       %J = StampCurrentSource(  J,  Ibr_eq, c,  b   );
       %J = StampCurrentSource(  J,  Icc_eq, e,  c   );
       
       Y = Stamp(   Y,  -gee,                   e,  e);
       Y = Stamp(   Y,  -gec,                   e,  c);
       Y = Stamp(   Y,  gee+gec,                e,  b);
       
       Y = Stamp(   Y,  -gce,                   c,  e);
       Y = Stamp(   Y,  -gcc,                   c,  c);
       Y = Stamp(   Y,  gce+gcc,                c,  b);
       
       Y = Stamp(   Y,  gee+gce,                b,  e);
       Y = Stamp(   Y,  gec+gcc,                b,  c);
       Y = Stamp(   Y,  -(gee+gec+gce+gcc),     b,  b);
       
       J = StampCurrentSource(  J,  -ie, e,  b   );
       J = StampCurrentSource(  J,  -ic, c,  b   );
       
       
       %Y = StampConductance(    Y,  gPiF,   b,  e   );
       %Y = StampConductance(    Y,  gPiR,   b,  c   );
       %Y = StampConductance(    Y,  go,     e,  c   );
       %Y = StampVCCS(           Y,  gmF,    b,  e,  c,  e 	);
       %Y = StampVCCS(           Y,  gmR,    b,  c,  e,  c   );
       
       
      %J = StampCurrentSource(  J,  Ibf_eq, e,  b   );
      %J = StampCurrentSource(  J,  Ibr_eq, c,  b   );
       %J = StampCurrentSource(  J,  Icc_eq, e,  c   );
       
     %  Y
     %  J
end%a

%% Fill the I matrix ##################################################


for j=1:numNode,
    for i=1:numD,
            
        n1 = Diodes(i).Node1 ;
        n2 = Diodes(i).Node2 ;
        if(n1~=0) v1=vNodeset(n1);
        else v1=0; end
        if(n2~=0) v2=vNodeset(n2);
        else v2=0; end
        
        %diode dc analysis working point       
        vd=v1-v2 ;          
        Id = Is * ( exp(vd/Vt)-1 ); %standard diode equation
        geq = Id/Vt;        % g_eq = i(v_n)/V_t         
        Ieq = Id - geq*vd ;
        
        if (Diodes(i).Node1==j),
            J(j)= J(j) - Ieq;
        elseif (Diodes(i).Node2==j),
            J(j)= J(j) + Ieq;
        end
    end
end

for j=1:numNode,
    for i=1:numI,
        if (Isource(i).Node1==j),
            if (symbolic==1)
            JS{j}=[JS{j} '-' Isource(i).Name];
            end
            J(j)= J(j) - Isource(i).Value;
        elseif (Isource(i).Node2==j),
            if (symbolic==1)
            JS{j}=[JS{j} '+' Isource(i).Name];
            end
            J(j)=J(j) + Isource(i).Value;
        end
    end
end

for i=1:numV,
    if (symbolic==1)
    JS{numNode+i} =  [JS{numNode+i} '+' Vsource(i).Name] ;
    end
    J(numNode+i) = J(numNode+i) + Vsource(i).Value ;
end

for i=1:numE,
    if (symbolic==1)
    JS{numNode+numV+i} = '0' ;
    end
    J(numNode+numV+i) = J(numNode+numV+i) + 0 ;
end
 %Linear Current-Controlled Current Sources
 %needs special consideration
for i=1:numF,
    %JS{numNode+numV+numE+i} = Fsource(i).Dependant ;
    %J(numNode+numV+numE+i) = J(numNode+numV+numE+i) + Vsource( Fsource(i).DependantID ).Value ;
    if (symbolic==1)
    JS{numNode+numV+numE+i} = [ JS{numNode+numV+numE+i} '+0' ] ;
    end
    J(numNode+numV+numE+i) = J(numNode+numV+numE+i) + 0 ; % Vsource( Fsource(i).DependantID ).Value ;
end
 
 %Linear Voltage-Controlled Current Sources
%for i=1:numG,
%    JS{numNode+numE+numF+i} = '0' ;
%    J(numNode+numE+numF+i) = 0 ;
%end
 
%Linear Current-Controlled Voltage Sources:
%needs special consideration
for i=1:numH,
    %JS{numNode+numV+numE+numF+i-1} = Hsource(i).Dependant ;
    %J(numNode+numV+numE+numF+i-1) = J(numNode+numV+numE+numF+i-1)+ Vsource( Hsource(i).DependantID ).Value ;
    if (symbolic==1)
    JS{numNode+numV+numE+numF+ i} = [JS{numNode+numV+numE+numF+ i} '+0'] ;
    end
    J(numNode+numV+numE+numF+ i) = J(numNode+numV+numE+numF+ i) + 0 ;
    
end
%The I matrix is done -----------------------------------------------


%% Fill the V matrix ##################################################
for i=1:numNode,
    V{i}=['V(' num2str(i) ')'];
end

for i=1:numV,
    V{numNode+i}=[sprintf('i(%d,%d)',  Vsource(i).Node1, Vsource(i).Node2)];
end


for i=1:numE,
    V{i+numNode+numV}=[sprintf('i(%d,%d)',  Esource(i).Node1, Esource(i).Node2)] ;
end
 %Linear Current-Controlled Current Sources
for i=1:numF,
    V{i+numNode+numE+numV} = [sprintf('i(%d,%d)',  Fsource(i).Node3, Fsource(i).Node4)];
end
 
 %Linear Voltage-Controlled Current Sources
%for i=1:numG,
%end
 %Linear Current-Controlled Voltage Sources:
for i=1:numH,
    %V{2*i-1+numNode+numE+numF+numV} =[sprintf('i(%d,%d)',  Hsource(i).Node3, Hsource(i).Node4)];
    V{i+numNode+numE+numF+numV} = [sprintf('i(%d,%d)',  Hsource(i).Node1, Hsource(i).Node2)];
    
     
end

%The V matrix is finished -------------------------------------------


%% Dependent sources
%Voltage-Controlled Voltage Sources
for i=1:numV,
    m=numNode+i-1;
    I=numNode+i;
    
    j=Vsource(i).Node1 ;
    jb=Vsource(i).Node2 ;
    
    if(j~=0),
        if (symbolic==1)
            YS{m+1, j} = [YS{m+1, j} '+1'];
            YS{j, I} = [YS{j, I} '+1'];
        end
        Y(m+1,j) = Y(m+1,j)+1;
        Y(j,I) = Y(j,I)+1;
    end
    
    if(jb~=0),
         if (symbolic==1)
            YS{m+1, jb} = [YS{m+1, jb} '-1'];
            YS{jb, I} = [YS{jb, I} '-1'];
         end
        
        Y(m+1,jb) = Y(m+1,jb)-1;
        Y(jb,I) = Y(jb,I)-1;
    end
    
end

for i=1:numE,
    m = numNode+numV + i -1 ;
    I= numNode+numV + i ;
    k=Esource(i).Node1 ;
    kb=Esource(i).Node2 ;
    j=Esource(i).Node3 ;
    jb=Esource(i).Node4 ;  
    u=Esource(i).Value ;
    if( j ~= 0 ),
        if (symbolic==1)
            YS{m+1,j} = [ YS{m+1,j} '-u' ];  
        end
        Y(m+1,j) = Y(m+1,j) - u ;
    end
    if( jb ~= 0 ),
        if (symbolic==1)
            YS{m+1,jb} = [ YS{m+1,jb} '+u'] ;
        end
        Y(m+1,jb) = Y(m+1,jb) + u ;
    end
    if( k ~= 0 ),
        if (symbolic==1)
        YS{m+1,k} = [ YS{m+1,k} '+1' ];
        YS{k, I} = [ YS{k, I} '+1'  ] ;
        end
        
        Y(m+1,k) = Y(m+1,k) +1 ;
        Y(k,I) = Y(k,I) +1 ;
    end
    if( kb ~= 0 ),
        if (symbolic==1)
        YS{m+1,kb} = [ YS{m+1,kb} '-1' ] ; 
        YS{kb, I} = [ YS{kb, I} '-1' ];
        end
        Y(m+1,kb) = Y(m+1,kb) - 1 ;
        Y(kb,I) = Y(kb,I) - 1 ;
    end
end
 %Linear Current-Controlled Current Sources
for i=1:numF,
    m = numNode+numV + numE + i -1 ;
    I= numNode+numV + numE +i ;
    ID = numNode + Fsource(i).DependantID ;
    k=Fsource(i).Node1 ;
    kb=Fsource(i).Node2 ;
    j=Fsource(i).Node3 ;
    jb=Fsource(i).Node4 ;  
    a=Fsource(i).Value ;
    %if( j ~= 0 ),
    %    YS{m+1,j} = [ YS{m+1,j} '+1' ];  
    %    Y(m+1,j) = Y(m+1,j) + 1 ;
        
    %    YS{j,I} = [ YS{j,I} '1' ];
    %    Y(j,I) = Y(j,I) +1 ; 
    %end
    %if( jb ~= 0 ),
    %    YS{m+1,jb} = [ YS{m+1,jb} '+1'] ;
    %    Y(m+1,jb) = Y(m+1,jb) -1 ;
    %    
    %    YS{jb,I} =  [ YS{jb,I} '-1' ];
    %    Y(jb,I) =  Y(jb,I)-1 ;
    %end
    
    %stamping the dependant voltage source
    if( k ~= 0 ),
        if (symbolic==1)
        YS{k, ID} = [ YS{k, ID} '+a'  ] ;
        end
        Y(k,ID) = Y(k,ID) +a ;
    end
    
    if( kb ~= 0 ),
        if (symbolic==1)
        YS{kb, ID} = [ YS{kb, ID} '-a' ];
        end
        Y(kb,ID) = Y(kb,ID) - a ;
    end
    
    %I1=a*I2:
    if (symbolic==1)
        YS{m+1, ID} = [ YS{m+1, ID} '-a'];
        YS{m+1, I} = [ YS{m+1, I} '1'];
    end
    
    Y(m+1,ID)=Y(m+1,ID)-a;
    Y(m+1,I)=Y(m+1,I)+1;
    
    
end
 
 %Linear Voltage-Controlled Current Sources
for i=1:numG,
    k=Gsource(i).Node1 ;
    kb=Gsource(i).Node2 ;
    j=Gsource(i).Node3 ;
    jb=Gsource(i).Node4 ;  
    g=Gsource(i).Value ;
    StampVCCS( Y, g, k, kb, j, jb );
end

 %Linear Current-Controlled Voltage Sources:
for i=1:numH,
    m = numNode+numV + numE + numF + i -1 ;
    I= numNode+numV + numE + numF + i  ;
    ID = numNode + Hsource(i).DependantID ;
    
    k=Hsource(i).Node1 ;
    kb=Hsource(i).Node2 ;
    j=Hsource(i).Node3 ;
    jb=Hsource(i).Node4 ;  
    r=Hsource(i).Value ;
    %if( j ~= 0 ),
    %    YS{m+1,j} = [ YS{m+1,j} '+1' ];  
    %    Y(m+1,j) = Y(m+1,j) +1 ;
        
    %    YS{j, I} = [ YS{j, I} '+1' ];  
    %    Y(j, I) = Y(j, I) +1 ;
    %end
    
    %if( jb ~= 0 ),
    %    YS{m+1,jb} = [ YS{m+1,jb} '-1'] ;
    %    Y(m+1,jb) = Y(m+1,jb) -1 ;
        
    %    YS{jb, I1} = [ YS{jb, I1} '-1' ];  
    %    Y(jb, I1) = Y(jb, I1) -1 ;
        
    %end
    if( k ~= 0 ),
        if (symbolic==1),
            YS{k, I} = [ YS{k, I} '+1'  ] ;
            YS{m+1,k} = [ YS{m+1,k} '+1' ];
        end
        
        Y(m+1,k) = Y(m+1,k) +1 ;
        
        
        Y(k,I) = Y(k,I) +1 ;
    end
    if( kb ~= 0 ),
        if (symbolic==1),
            YS{m+1,kb} = [ YS{m+1,kb} '-1' ] ; 
            YS{kb, I} = [ YS{kb, I1} '-1' ];
        end
        
        Y(m+1,kb) = Y(m+1,kb) - 1 ;
        Y(kb,I) = Y(kb,I1) - 1 ;
    end
    
    if( ID~=0),
        if (symbolic==1),
            YS{m+1, ID} = [YS{m+1,I} '-r'];
        end
        Y(m+1, ID) = Y(m+1,I) -r ;
    end
    
end




%% DC-Analysis
if (numD>0),
    disp('DC iterations ...');
end 



Y
J
%first approach:
dcPointVoltage = Y\J ; % this is faster than dcPointVoltage = inv(Y)*J ;

%   [L,U] = lu(Y);
%   y = L\J;
%   x = U\y;    %internally, matlab uses LU factorization in Y\J.
%Second Approach: Using LU Factorization and Guassian Elimination


%% Printing out some data:
    if( diff == 1 ),
        [size t] = size(dcPointVoltage);
        for i=1:size
            disp( sprintf('[DC] %s = ', V{i} ) )
            dcPointVoltage(i)
        end
    end
    
    dcPointVoltage
    d = vNodeset - dcPointVoltage(1:numNode) ;
    vNodeset = dcPointVoltage(1:numNode) ;
    diff = sqrt(d'*d)
    iteration = iteration + 1; 
    disp( sprintf('iteration %d, difference=%f', iteration, diff))
end

%invA=inv(A);
 
%solve and store results
%for i=2:length(t),
%    b=[0 0 Vin(i)]';
%    x=invA*b;
%    V1(i)=x(1);  V2(i)=x(2);  V3(i)=x(3);
%end


disp(sprintf('[info] Elapsed time = %g seconds.\n',toc));
beep;