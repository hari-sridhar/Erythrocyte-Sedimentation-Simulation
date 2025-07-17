function [v, So, L, Adj, source_node_interparticle, node_currents] = computeNodeVoltages(nodes, edges, Res, N, source_nodes, sink_nodes, num_edges, num_original_nodes)

% Source and Voltage Calculation

Adj = zeros(num_original_nodes, num_original_nodes);
for k = 1:num_edges
    i = edges(k, 1); 
    j = edges(k, 2);
    if nodes(j, 2) > nodes(i, 2)
        Adj(i, j) = 1/Res(i, j);
    end
end

D = diag(sum(Adj, 2));
L = D - Adj;


So = zeros(num_original_nodes, 1);
source_node_interparticle = false(num_original_nodes, 1);
node_currents = zeros(num_original_nodes, 1);

for i = 1:N
    node_1 = 4*(i-1) + 1;
    has_interparticle = any((edges(:,1) == node_1 | edges(:,2) == node_1) & ...
                                ((ceil(edges(:,1)/4) ~= i) | (ceil(edges(:,2)/4) ~= i)));        
    if ~has_interparticle
        So(node_1) = 1;
        node_currents(node_1) = 1;
    else
        source_node_interparticle(node_1) = true;
    end
end

So(source_nodes) = 1;
So(sink_nodes) = -1;  % You are correct to make all sinks -1

% Solve
v = pinv(L)* So;