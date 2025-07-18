function edges = generateParticleEdges(N)
% Intraparticle Edges

edges = [];
for i = 1:N
    edges = [edges; 
        4*(i-1)+1, 4*(i-1)+2;  % node1-node2
        4*(i-1)+1, 4*(i-1)+3;  % node1-node3
        4*(i-1)+2, 4*(i-1)+4;  % node2-node4
        4*(i-1)+3, 4*(i-1)+4]; % node3-node4
end
