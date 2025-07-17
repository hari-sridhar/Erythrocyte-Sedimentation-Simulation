function animateParticleStreamline(X_history, Y_history, R, colors, ...
                                    x_min, x_max, nodes_all, edge_currents_all, alpha)

% Padding and axis limits
padding = 10e-6;
y_minimum_simulate = min(Y_history(:));
y_maximum_simulate = max(Y_history(:));
[N, num_time_steps] = size(X_history);
figure;

% Define colormap
cmap = jet(256);

% ✅ Robust max current calculation
max_I_global = 0;
for t = 1:length(edge_currents_all)
    E = edge_currents_all{t};
    if isempty(E) || size(E, 2) < 3
        continue;
    end
    local_max = max(abs(E(:,3)));
    if local_max > max_I_global
        max_I_global = local_max;
    end
end
if max_I_global == 0
    max_I_global = 1e-12;  % Fallback to prevent divide by zero
end

% Animation loop
for t = 1:num_time_steps
    clf;
    hold on;
    axis equal;

    % Background plasma
    fill([x_min x_max x_max x_min], [y_minimum_simulate y_minimum_simulate y_maximum_simulate  y_maximum_simulate ], ...
         [1 0.85 0.85], 'EdgeColor', 'none');

    % Particles
    for j = 1:N
        rectangle('Position', [X_history(j, t) - R(j), Y_history(j, t) - R(j), 2 * R(j), 2 * R(j)], ...
                  'Curvature', [1, 1], 'EdgeColor', colors(j, :), ...
                  'FaceColor', 'w', 'LineWidth', 1.5);
        text(X_history(j, t), Y_history(j, t), num2str(j), ...
             'FontSize', 8, 'Color', 'k', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    % Axis limits (zoom to particle region with padding)
    xmin_frame = min(X_history(:, t)) - padding;
    xmax_frame = max(X_history(:, t)) + padding;
    ymin_frame = min(Y_history(:, t)) - padding;
    ymax_frame = max(Y_history(:, t)) + padding;
    xlim([xmin_frame, xmax_frame]);
    ylim([ymin_frame, ymax_frame]);

    % Edge currents and arrows
    nodes = nodes_all{t};
    edge_currents = edge_currents_all{t};

    if ~isempty(edge_currents) && size(edge_currents, 2) >= 3
        for k = 1:size(edge_currents, 1)
            i = edge_currents(k, 1);
            j = edge_currents(k, 2);
            I = edge_currents(k, 3);

            xi = nodes(i, 1); yi = nodes(i, 2);
            xj = nodes(j, 1); yj = nodes(j, 2);

            % Flip direction if I < 0
            if I < 0
                [xi, xj] = deal(xj, xi);
                [yi, yj] = deal(yj, yi);
            end

            % Recompute direction and midpoint
            dx = xj - xi;
            dy = yj - yi;
            if dx == 0 && dy == 0
                continue;
            end
            xm = (xi + xj) / 2;
            ym = (yi + yj) / 2;

            % Rotation angle
            angle_deg = atan2d(dy, dx);
            if angle_deg > 90
                angle_deg = angle_deg - 180;
            elseif angle_deg < -90
                angle_deg = angle_deg + 180;
            end

            % Color by current magnitude
            idx = round(255 * min(abs(I) / max_I_global, 1)) + 1;
            arrow_color = cmap(idx, :);

            % Draw arrow
            scale = 0.2;
            quiver(xi, yi, dx * scale, dy * scale, 0, ...
                   'Color', arrow_color, 'LineWidth', 1.5, 'MaxHeadSize', 2);

            % Draw rotated annotation at edge midpoint
            text(xm, ym, sprintf('%.2e', abs(I)), ...
                'FontSize', 8, 'Color', 'k', ...
                'HorizontalAlignment', 'center', ...
                'Rotation', angle_deg);
        end
    end

    % Walls
    xline(x_min, '--k', 'LineWidth', 2);
    xline(x_max, '--k', 'LineWidth', 2);

    % Labels and axis
    title(sprintf('Edge current Flow with Particle Motion (with Tunable Parameter set at %.2e)', alpha) );
    xlabel('X-axis'); ylabel('Y-axis');

    drawnow;
    pause(0.01);
end
