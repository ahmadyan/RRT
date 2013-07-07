function core = init(mode, probe1, probe2, probe3)
format long e
core.numElem=0;  %Number of passive elements.
core.numV=0;     %Number of independent voltage sources
core.numO=0;     %Number of op amps
core.numI=0;     %Number of independent current sources
core.numD=0;

core.numE=0;
core.numF=0;
core.numG=0;
core.numH=0;
core.numQ=0;

core.numNode=0;  %Number of nodes, not including ground (node 0).
core.node1 =0;
core.node2 =0;

core.components = java.util.Hashtable ;
core.nodes = java.util.Hashtable ;
core.components.clear;
core.nodes.clear;
core.nodes.put ('0', 0); %Adding ground
core.numNode = 1 ;
core.operatingPoint=[];

core.isNonlinear=0;

%decleras the operation mode
% DC Analysis=0
% AC Analysis=1
% Transient Analysis=2

core.mode=mode;    

%ac analysis
core.freq=0;
core.fmin=10e3;
core.fmax=100e9;
core.step=1000;

%transient simulation parameters
core.t0 = 0;
core.tf = 200e-6;
core.dt = 0.1E-6;

core.dcIteration=0;

%to find probe#, write core.nodes. Then add a node number here.
%for example: {input2=5.0, b=4.0, input1=6.0, cbias=2.0, out=3.0, 0=0.0, supply=1.0}
%so adding a core.probe=3 will monitor the out pin.
%probes for amp_tran

core.probe=probe1;
core.probe2=probe2;
core.probe3=probe3;


core.enableDebug=0;

end