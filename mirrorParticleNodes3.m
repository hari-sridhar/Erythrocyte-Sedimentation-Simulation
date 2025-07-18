function [nodes_ext, labels_ext, num_nodes_ext, mirrored_particles] = mirrorParticleNodes(num_original_nodes, x_min, x_max, nodes)


% Original labels
labels_orig = (1:num_original_nodes)';

% Mirror across x_min
nodes_left = nodes;
nodes_left(:,1) = 2*x_min - nodes(:,1);
labels_left = -labels_orig;

% Mirror across x_max
nodes_right = nodes;
nodes_right(:,1) = 2*x_max - nodes(:,1);
labels_right = labels_orig + num_original_nodes;

% Combine nodes and labels
nodes_ext  = [nodes;      nodes_left;      nodes_right];
labels_ext = [labels_orig; labels_left;     labels_right];

% Output mirrored positions for diagnostic if needed
t.mirrored_particles.left  = nodes_left;
t.mirrored_particles.right = nodes_right;
mirrored_particles = t;

% Update count
num_nodes_ext = size(nodes_ext,1);
end
