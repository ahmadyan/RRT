function [ nodes ] = rrt2mati( inputFileName )
fid = fopen( inputFileName, 'r');

fgets(fid); %rrt
line = fgets(fid); %number of dimension
dim = sscanf(line,'%d');

line = fgets(fid); %number of variations
var = sscanf(line,'%d');

line = fgets(fid); %number of samples
samples = sscanf(line,'%d')-1;

min = zeros(1, dim);
max = zeros(1, dim);
for i=1:dim,
    line = fgets(fid); %minimum value of ith dim
    min(i) = sscanf(line,'%f');
    line = fgets(fid); %minimum value of ith dim
    max(i) = sscanf(line,'%f');
end

nodes = zeros(dim+var+2, samples+1);
%while ~feof(fid)
for i=0:samples,
    line = fgets(fid);
    %1 0 0.647759 -0.00207225 5e-009 
    nodeRawData = sscanf(line, '%f');
    nodes(:, i+1) = nodeRawData;
    nodes(1, i+1) = nodes(1, i+1)+1;    %in matlab, arrays are indexed at 1, instead of 0 (c++)
    nodes(2, i+1) = nodes(2, i+1)+1;
end

end
