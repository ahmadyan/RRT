function golden=getDev(data, k, dim)
j=2;
for i=1:size(data,2)
    if(mod(i-1, k)==1),
        golden(:,j) = data(1:dim, i);
        j=j+1;
    end
end