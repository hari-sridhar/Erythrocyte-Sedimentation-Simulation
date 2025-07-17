function plotParticleDynamics(X_history, Y_history, Vx_history, Vy_history, time_steps, colors, N, alpha)

% --- Plotting the trajectory ---
figure('Position', [100, 100, 900, 700]);
hold on;

% Build strings for display below legend
info_lines = strings(N, 1);  % Preallocate array of strings

for j = 1:N
    plot(X_history(j, :), Y_history(j, :), 'LineWidth', 2, ...
        'DisplayName', sprintf('Particle %d', j), ...
        'Color', colors(j, :));

    % Construct info string
    info_lines(j) = sprintf('P%d: Start(%.2e, %.2e), End(%.2e, %.2e)', ...
        j, X_history(j, 1), Y_history(j, 1), X_history(j, end), Y_history(j, end));
end

legend();
xlabel('X Position');
ylabel('Y Position');
title(sprintf('Trajectories of Particles (Phase 3 with Tunable Parameter set at %.2e)', alpha));
grid on;
axis equal;

% Add annotation text box below plot
annotation_text = strjoin(info_lines, '\n');  % Join all lines with newline
dim = [0.15, 0.01, 0.8, 0.1];  % [x y w h] relative to figure
annotation('textbox', dim, 'String', sprintf(annotation_text), ...
           'FitBoxToText', 'on', 'EdgeColor', 'none', ...
           'FontSize', 8, 'Interpreter', 'none');

% --- X-Velocity Plot ---
figure;
hold on;
for j = 1:N

    plot(time_steps, Vx_history(j, :), 'LineWidth', 2, ...
         'DisplayName', sprintf('Particle %d', j), ...
         'Color', colors(j, :));

        % Construct info string
    info_lines(j) = sprintf('P%d: Start(%.2e), End(%.2e)', ...
        j, Vx_history(j, 1), Vx_history(j, end));

end
legend();
xlabel('Time');
ylabel('X Velocity');
title(sprintf('X-Velocity of Particles (Phase 3 with Tunable Parameter set at %.2e)', alpha));
grid on;

% Add annotation text box below plot
annotation_text = strjoin(info_lines, '\n');  % Join all lines with newline
dim = [0.15, 0.01, 0.8, 0.1];  % [x y w h] relative to figure
annotation('textbox', dim, 'String', sprintf(annotation_text), ...
           'FitBoxToText', 'on', 'EdgeColor', 'none', ...
           'FontSize', 8, 'Interpreter', 'none');

% --- Y-Velocity Plot ---
figure;
hold on;
for j = 1:N
    plot(time_steps, Vy_history(j, :), 'LineWidth', 2, ...
         'DisplayName', sprintf('Particle %d', j), ...
         'Color', colors(j, :));

            % Construct info string
    info_lines(j) = sprintf('P%d: Start(%.2e), End(%.2e)', ...
        j, Vy_history(j, 1), Vy_history(j, end));


end
legend();
xlabel('Time');
ylabel('Y Velocity');
title(sprintf('Y-Velocity of Particles (Phase 3 with Tunable Parameter set at %.2e)', alpha));
grid on;

% Add annotation text box below plot
annotation_text = strjoin(info_lines, '\n');  % Join all lines with newline
dim = [0.15, 0.01, 0.8, 0.1];  % [x y w h] relative to figure
annotation('textbox', dim, 'String', sprintf(annotation_text), ...
           'FitBoxToText', 'on', 'EdgeColor', 'none', ...
           'FontSize', 8, 'Interpreter', 'none');

hold off;
