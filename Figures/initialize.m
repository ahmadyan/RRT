function [Golden,deviation,upper,lower]=initialize(data,goldenvalue)


 Golden=getGolden(goldenvalue);
 Deviation=getDev(Golden,data);
 [upper,lower]=findEnv(Deviation);
 
 
 
 

