function [nodes, node_labels, num_original_nodes, highest_particle, particle_voronoi_neighbors] = generateParticleNodes(X, Y, R, N)

% Delaunay Triangualtion, Voronoi Neighbours, Nodes and Edges
particle_centers = [X,Y];
[~, highest_particle] = max(Y);


    dt_particles = delaunayTriangulation(particle_centers(:,1), particle_centers(:,2));
    tri = dt_particles.ConnectivityList;
    particle_voronoi_neighbors = zeros(N, N);
    for k = 1:size(tri, 1)
        for a = 1:3
            for b = a+1:3
                i = tri(k, a);
                j = tri(k, b);
                particle_voronoi_neighbors(i, j) = 1;
                particle_voronoi_neighbors(j, i) = 1;
            end
        end
    end

get_nodes = @(x, y, r) [x, y - r; x - r, y; x + r, y; x, y + r];
nodes = []; node_labels = [];
for i = 1:N
    new_nodes = get_nodes(particle_centers(i,1), particle_centers(i,2), R(i));
    nodes = [nodes; new_nodes];
    node_labels = [node_labels; 4*(i-1)+1; 4*(i-1)+2; 4*(i-1)+3; 4*(i-1)+4];
    
end

% Save original number of nodes
num_original_nodes = size(nodes,1);
