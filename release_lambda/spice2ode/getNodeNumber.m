function nodeID = getNodeNumber( nodes, nodeName )
%GETRNODENUMBER get NodeID for a given node name
%   nodes is a hashtable (from java.util.hashtable that contains the
%   corresponding node id for each node name.
    if( nodes.containsKey( nodeName ) == 1 )
        nodeID = nodes.get( nodeName ) ;
    else
        nodes.put( nodeName, size(nodes) ) ;
        nodeID = nodes.get( nodeName ) ;
    end
end

