function [highest_4th_node, sink_nodes, source_nodes, num_edges] = generateSource(highest_particle, N, edges)
% Sink and Source

highest_4th_node = 4 * highest_particle;

sink_nodes = [];
source_nodes = [];
%particle_connections = zeros(N, N);  % Track if particle i connected to j


for i = 1:N
    source = 4 * (i-1)+1;  % 4th node of particle i
    sink = 4*i;
    sink_nodes = [sink_nodes; sink];   
    source_nodes = [source_nodes; source];
end

num_edges = size(edges,1);