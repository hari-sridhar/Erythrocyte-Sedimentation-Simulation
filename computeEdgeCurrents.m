function [incoming_current, edge_currents, I_values, total_current] = computeEdgeCurrents (edges_ascending, num_original_nodes, Res, v, source_node_interparticle, node_currents, num_edges)

% Current calculation

incoming_current = zeros(num_original_nodes, 1);
edge_currents = zeros(num_edges, 3);
I_values = zeros(num_edges, 3);  

for k = 1:num_edges
    i = edges_ascending(k, 1); j = edges_ascending(k, 2);
    I = (v(i) - v(j)) / Res(i, j);
    I_values(k,:) = [i, j, I];

    if source_node_interparticle(i)
        incoming_current(j) = (incoming_current(j) + incoming_current(i) + I)/2;
        edge_currents(k,:) = [i, j, (incoming_current(i) + I) / 2];
    else
        incoming_current(j) = incoming_current(j) + incoming_current(i) + I;
        edge_currents(k,:) = [i, j, incoming_current(i) + I];
    end
end

% incoming_current = zeros(num_nodes, 1);
% edge_currents = zeros(num_edges, 3);
% I_values = zeros(num_edges, 3);  
% 
% for k = 1:num_edges
%     i = edges_ascending(k, 1);
%     j = edges_ascending(k, 2);
% 
%     % Calculate edge current
%     I = (v(i) - v(j)) / R(i, j);
%     I_values(k,:) = [i, j, I];
% 
%     % Direct current addition
%     edge_currents(k,:) = [i, j, incoming_current(i) + I];
% 
%     % Update incoming current at j
%     incoming_current(j) = incoming_current(j) + incoming_current(i) + I;
% end

total_current = node_currents + incoming_current;