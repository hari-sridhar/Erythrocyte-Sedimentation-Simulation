clc; clear; close all;

%% User Input

% Define the matrix size
N = input("Enter the number of particles: ");

% Blood Parameters
eta = 0.0012; % Viscosity of Blood Plasma
m = 27e-12; % Mass of Particles
rho = 1025; % Density of Plasma
g = 9.81; % Gravitational constant
alpha = 1e-8; %Tuning Parameter


% Time parameters
t_start = 0;
t_end = 0.25;  % Total simulation time
dt = 0.001;    % Time step
time_steps = t_start:dt:t_end;
num_time_steps = length(time_steps);



% Initialize position history matrices for N particles over time
X_history = zeros(N, num_time_steps);
Y_history = zeros(N, num_time_steps);

% Define the domain Boundaries 
x_min = 0e-6; x_max = 150e-6;
y_min = 0e-6; y_max = 150e-6;


% Initial Positioning
if N == 5
    X = [7e-6; 50e-6; 70e-6; 30e-6; 62e-6];
    Y = [110e-6; 110e-6; 110e-6; 110e-6; 110e-6];
    R = [4.07e-06; 3.84e-06; 4.00e-06; 3.56e-6; 3.78e-6];

elseif N==3
    X = [7e-6; 40e-6; 25e-6; ];
    Y = [110e-6; 110e-6; 110e-6 ];
    R = [4.07e-06; 3.84e-06; 4.00e-06];

elseif N==2
    X = [15e-6; 40e-6];
    Y = [110e-6; 110e-6];
    R = [3.6e-6; 3.8e-6];
else
    X = (1.05*x_min) + (0.95*(x_max - x_min)) * rand(N,1);
    
    Y  = (0.95*y_min) + (0.95*(y_min - y_max))* zeros(N,1);
    
    % Create a Radius matrix
    R = 3.6e-6 + (4.1e-6 - 3.6e-6) * rand(N, 1); % Random radii for particles

    separation = 2.1 * max(R);  % Minimum center-to-center separation

    for k = 2:N
    valid = false;
        while ~valid
        X(k) = x_min + (x_max - x_min)*rand;
        % check horizontal spacing only:
        valid = all(abs(X(k) - X(1:k-1)) >= separation);
        end
    end
end


% Initialize velocity history matrices for N particles over time
Vx_history = zeros(N, num_time_steps);
Vy_history = zeros(N, num_time_steps);

% Initialize the unit vector matrices for x and y components
ex = zeros(N, N);  % Matrix to hold x-components of unit vectors
ey = zeros(N, N);  % Matrix to hold y-components of unit vectors

D = zeros(N);

% Assign initial positions to the first column of X_history and Y_history
X_history(:, 1) = X;
Y_history(:, 1) = Y;

%Create a Kd matrix (Drag coefficient) 
KD = 6*pi*eta*R;

% Define volume
vol = (4/3) * pi * R.^3;

max_saves   = ceil(num_time_steps/10);
Res_history = cell(max_saves,1);
save_idx    = 1;

% Define random colors for each particle
colors = rand(N, 3);

% Set figure properties
%figure('Position', [100, 100, 1540, 1000]); % Large figure size



%% Simulation Loop (No Domain Boundaries)

for t = 2:num_time_steps


    %% Node Computation

%% Delaunay Triangualtion, Voronoi Neighbours, Nodes and Edges
[nodes, node_labels, num_original_nodes, highest_particle, particle_voronoi_neighbors] = generateParticleNodes(X, Y, R, N);


%% Mirror nodes across walls
%[node_labels, num_nodes, mirrored_particles] = mirrorParticleNodes(num_original_nodes, N, x_min, x_max, nodes);

[nodes_ext, labels_ext, num_nodes_ext, mirrored_particles] = mirrorParticleNodes3(num_original_nodes, x_min, x_max, nodes);

edges = generateParticleEdges(N);

%% Sink Nodes and Edges

[highest_4th_node, sink_nodes, source_nodes, num_edges] = generateSource(highest_particle, N, edges);

%% Prepare edge list such that i is lower y than j
edges_ascending = getEdgesAscending(edges,nodes);

%% Update particle interactions
    for p = 1:N
        for q = p+1:N  % Only calculate for p < q to avoid redundant calculations
            
            % Distance between particles
            dx = X_history(q, t-1) - X_history(p, t-1);
            dy = Y_history(q, t-1) - Y_history(p, t-1);
            D(p, q) = sqrt(dx^2 + dy^2) - (R(p) + R(q)); 
            D(q, p) = D(p, q);
            
            % Calculate unit vector components if the distance is non-zero
            if D(p, q) ~= 0
                ex(p, q) = dx / D(p, q);
                ey(p, q) = dy / D(p, q);
                ex(q, p) = -ex(p, q);  % Symmetry
                ey(q, p) = -ey(p, q);
            end
        end
    end

    % Calculate P, Q, S matrices
    P = zeros(N); Q = zeros(N); S = zeros(N);
    for j = 1:N
        for k = 1:N
            if j ~= k
                P(j, k) = (6 * pi * eta * (R(j) + R(k))^2 / 4) * ex(j, k)^2 / D(j, k);
                P(k, j) = P(j, k);
                Q(j, k) = (6 * pi * eta * (R(j) + R(k))^2 / 4) * ex(j, k) * ey(j, k) / D(j, k);
                Q(k, j) = Q(j, k);
                S(j, k) = (6 * pi * eta * (R(j) + R(k))^2 / 4) * ey(j, k)^2 / D(j, k);
                S(k, j) = S(j, k);
            end
        end
    end

    %disp("Size of nodes matrix:"); disp(size(nodes));
    %disp("Maximum index in edges matrix:"); disp(max(edges(:)));
    %disp(edges)


%% Triangulation

%[Res, triangle_info] = computeResistanceMatrix(nodes, node_labels, num_original_nodes, edges, particle_voronoi_neighbors, R, eta, y_max);

%[triangle_info, Res] = compute_resist_triangle(nodes, node_labels, N, eta, y_max);

[triangle_info, Res] = compute_wall_edge_resistance(nodes_ext, labels_ext, nodes, node_labels, N, eta, y_max);

% 


%% Voltage calculation

[v, So, L, Adj, source_node_interparticle, node_currents] = computeNodeVoltages(nodes, edges, Res, N, source_nodes, sink_nodes, num_edges, num_original_nodes);

%% Current calculation
[incoming_current, edge_currents, I_values, total_current] = computeEdgeCurrents(edges_ascending, num_original_nodes , Res, v, source_node_interparticle, node_currents, num_edges);

% Update node positions based on particle positions
    nodes = zeros(N*4, 2);
    for j = 1:N
        x = X_history(j, t-1);
        y = Y_history(j, t-1);
        rj = R(j);
        base = 4*(j-1);
        nodes(base+1, :) = [x, y - rj];
        nodes(base+2, :) = [x - rj, y];
        nodes(base+3, :) = [x + rj, y];
        nodes(base+4, :) = [x, y + rj];
    end
    
nodes_all{t} = nodes;
edge_currents_all{t} = edge_currents(:,:);

%% KD Calculation

f_para = wall_lubrication_correction(X, R, N, x_min, x_max); 

KD    = 6*pi*eta .* R .* f_para;   % or however you combine KD & f_para
for j = 1:N
    nodes_j   = 4*(j-1) + (1:4);
    I_left    = incoming_current(nodes_j(2));
    I_right   = incoming_current(nodes_j(3));
    delta_I   = I_left - I_right;

    P(j,j) = -( KD(j) + alpha*delta_I + sum( P(j,:) ) );
    Q(j,j) = -(    sum( Q(j,:) ) );
    S(j,j) = -( KD(j) + sum( S(j,:) ) );
end

    % Solve for velocities
    A = [P, Q; Q, S];
    V = (m - rho * vol) * g;
    B = [zeros(N, 1); V];
    velocities = linsolve(A, B); 

    % Extract velocities
    Vx = velocities(1:N);
    Vy = velocities(N+1:end);

    % Update positions based on velocities
    X_new = X_history(:, t-1) + Vx * dt;
    Y_new = Y_history(:, t-1) + Vy * dt;

    % Store updated positions and velocities
    X_history(:, t) = X_new;
    Y_history(:, t) = Y_new;
    Vx_history(:, t) = Vx;
    Vy_history(:, t) = Vy;
end

% Graph setup
G = digraph(edge_currents(:,1), edge_currents(:,2), edge_currents(:,3));
x = nodes(1:num_original_nodes, 1);
y = nodes(1:num_original_nodes, 2);

fprintf("---- Res Element-wise View ----\n");

    for i = 1:size(Res,1)
        for j = i+1:size(Res,2)
            if Res(i,j) ~= inf && i ~= j && mod(j,4) ~= 1
            fprintf('Res(%d,%d) = %.4e\n', i, j, Res(i,j));
            end
        end
    end

%% Animation (Particles Moving Freely)

%animateParticles(X_history, Y_history, R, colors, x_min, x_max)

animateParticleStreamline(X_history, Y_history, R, colors, ...
                                    x_min, x_max, nodes_all, edge_currents_all, alpha)

%% Plotting the trajectory

plotParticleDynamics(X_history, Y_history, Vx_history, Vy_history, time_steps, colors, N, alpha)      

save('phase3_data.mat', 'X_history', 'Y_history');
