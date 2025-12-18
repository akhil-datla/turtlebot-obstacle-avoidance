% simulate_hybrid_two_turtlebots
%   Virtual simulation of two point-mass "Turtlebots" using the hybrid
%   Lyapunov-based obstacle avoidance controller from Sanfelice et al. (ACC
%   2006). This script implements the flow/jump sets with the
%   synergistic potentials (Eq. 11) and the distance functions (Eq. 10).

clear; clc; close all;

%% Geometry and potential parameters (meters)
% Per Example 5.2 in Sanfelice et al. (ACC 2006):
% - Single target at (3, 0)
% - Obstacle at (1, 0) with radius 1/(20*sqrt(2))
% - Mode q determines PATH around obstacle, not destination

% Single target (both robots go here via different paths)
% Per Example 5.2 in Sanfelice et al. (ACC 2006): target at (3, 0)
params.xt = 3.0;
params.yt = 0.0;

% Obstacle parameters (matching paper)
params.obs_cx = 1.0;
params.obs_cy = 0.0;
params.obs_R  = 1/(20*sqrt(2));  % δ = 1/(20*sqrt(2)) ≈ 0.0354 from paper

% Wedge angle for mode regions (per Figure 5)
params.wedge_angle = pi/4;  % 45 degrees

% Hybrid controller parameters
hyb.mu     = 1.2;   % hysteresis factor (>1)
hyb.lambda = 0.05;  % noise margin
hyb.gamma  = 0.01;  % small relaxation near the goal

dt       = 0.01;
T_final  = 30.0;
noise_std = 0.0;  % set >0 to study robustness (controller always uses y)

%% Initial conditions for two virtual robots
% Both robots navigate to the SAME target at (3, 0)
% Mode determines which PATH around the obstacle (upper vs lower)

% Robot 1: starts above, takes upper path (mode 1)
x1 = [0; 0.5];
q1 = 1;

% Robot 2: starts below, takes lower path (mode 2)
x2 = [0; -0.5];
q2 = 2;

%% Preallocate logging arrays
N = floor(T_final / dt);
t      = zeros(1, N);
x1_log = zeros(2, N);
x2_log = zeros(2, N);
q1_log = zeros(1, N);
q2_log = zeros(1, N);
V_active_log_1 = zeros(2, N); % potentials per robot
V_active_log_2 = zeros(2, N);
jump_log = zeros(2, N);
dist1_log = zeros(1, N); % distance to goal for robot 1
dist2_log = zeros(1, N); % distance to goal for robot 2

% Convergence parameters
goal_tolerance = 0.05; % meters
converged_1 = false;
converged_2 = false;

%% Simulation loop
for k = 1:N
    t(k) = (k-1) * dt;

    [x1_next, q1_next, info1] = hybrid_controller_step(x1, q1, params, hyb, dt, noise_std);
    [x2_next, q2_next, info2] = hybrid_controller_step(x2, q2, params, hyb, dt, noise_std);

    x1 = x1_next; q1 = q1_next;
    x2 = x2_next; q2 = q2_next;

    % Log trajectories and diagnostics
    x1_log(:, k) = x1;
    x2_log(:, k) = x2;
    q1_log(k) = q1;
    q2_log(k) = q2;
    V_active_log_1(:, k) = info1.V_y;
    V_active_log_2(:, k) = info2.V_y;
    jump_log(:, k) = [info1.jumped; info2.jumped];
    
    % Compute distance to goal (single target)
    goal = [params.xt; params.yt];
    dist1_log(k) = norm(x1 - goal);
    dist2_log(k) = norm(x2 - goal);
    
    % Check convergence
    if ~converged_1 && dist1_log(k) < goal_tolerance
        converged_1 = true;
        fprintf('Robot 1 reached goal at t = %.2f s\n', t(k));
    end
    if ~converged_2 && dist2_log(k) < goal_tolerance
        converged_2 = true;
        fprintf('Robot 2 reached goal at t = %.2f s\n', t(k));
    end
    
    % Early termination if both converged
    if converged_1 && converged_2
        fprintf('Both robots converged. Stopping simulation.\n');
        % Trim arrays
        t = t(1:k);
        x1_log = x1_log(:, 1:k);
        x2_log = x2_log(:, 1:k);
        q1_log = q1_log(1:k);
        q2_log = q2_log(1:k);
        V_active_log_1 = V_active_log_1(:, 1:k);
        V_active_log_2 = V_active_log_2(:, 1:k);
        dist1_log = dist1_log(1:k);
        dist2_log = dist2_log(1:k);
        break;
    end
end

%% Plot trajectories with obstacle and goal
figure; hold on; axis equal;
th = linspace(0, 2*pi, 200);
plot(params.obs_cx + params.obs_R*cos(th), params.obs_cy + params.obs_R*sin(th), 'k-', 'LineWidth', 2);
plot(params.xt, params.yt, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(x1_log(1, :), x1_log(2, :), 'b-', 'LineWidth', 1.5);
plot(x2_log(1, :), x2_log(2, :), 'g-', 'LineWidth', 1.5);
legend('Obstacle', 'Goal', 'Robot 1 (mode 1 start)', 'Robot 2 (mode 2 start)');
xlabel('x [m]'); ylabel('y [m]');
title('Hybrid obstacle avoidance (virtual Turtlebots)');

%% Plot mode evolution showing jump logic (hysteresis)
figure;
subplot(2,1,1);
stairs(t, q1_log, 'b', 'LineWidth', 1.5);
ylabel('q_1'); grid on;
title('Mode evolution with hysteresis (\mu for flow, \mu-\lambda for jumps)');
subplot(2,1,2);
stairs(t, q2_log, 'g', 'LineWidth', 1.5);
ylabel('q_2'); xlabel('Time [s]'); grid on;

%% Plot potentials versus time at the measurement points
figure;
subplot(2,1,1);
plot(t, V_active_log_1(1, :), 'b-', 'LineWidth', 1.2); hold on;
plot(t, V_active_log_1(2, :), 'r--', 'LineWidth', 1.2);
legend('V_1(y_1)', 'V_2(y_1)'); ylabel('Robot 1 potentials'); grid on;
subplot(2,1,2);
plot(t, V_active_log_2(1, :), 'g-', 'LineWidth', 1.2); hold on;
plot(t, V_active_log_2(2, :), 'm--', 'LineWidth', 1.2);
legend('V_1(y_2)', 'V_2(y_2)'); ylabel('Robot 2 potentials'); xlabel('Time [s]'); grid on;

%% Plot distance to goal versus time
figure;
subplot(2,1,1);
plot(t, dist1_log, 'b-', 'LineWidth', 1.5); hold on;
plot([t(1), t(end)], [goal_tolerance, goal_tolerance], 'r--', 'LineWidth', 1);
ylabel('Distance [m]'); grid on; title('Robot 1: Distance to goal');
legend('Distance', 'Goal tolerance');
subplot(2,1,2);
plot(t, dist2_log, 'g-', 'LineWidth', 1.5); hold on;
plot([t(1), t(end)], [goal_tolerance, goal_tolerance], 'r--', 'LineWidth', 1);
ylabel('Distance [m]'); xlabel('Time [s]'); grid on; title('Robot 2: Distance to goal');
legend('Distance', 'Goal tolerance');

% Print final statistics
fprintf('\n=== Simulation Summary ===\n');
fprintf('Final distance to goal:\n');
converged_1_str = 'false'; if converged_1, converged_1_str = 'true'; end
converged_2_str = 'false'; if converged_2, converged_2_str = 'true'; end
fprintf('  Robot 1: %.4f m (converged: %s)\n', dist1_log(end), converged_1_str);
fprintf('  Robot 2: %.4f m (converged: %s)\n', dist2_log(end), converged_2_str);
fprintf('Total mode jumps:\n');
fprintf('  Robot 1: %d\n', sum(jump_log(1, :)));
fprintf('  Robot 2: %d\n', sum(jump_log(2, :)));
fprintf('========================\n');

disp('Simulation complete.');
