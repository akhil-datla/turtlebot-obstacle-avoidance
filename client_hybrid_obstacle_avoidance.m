% client_hybrid_obstacle_avoidance
%   Real-time hybrid obstacle avoidance for two Turtlebots using the
%   Lyapunov-based controller from Sanfelice et al. (ACC 2006).
%   The controller always uses the measured position y (Vicon) to compute
%   V_q(y), d_q(y), and u = -\nabla V_q(y). The vector u is mapped to unicycle
%   commands via heading alignment every cycle.

clear; clc; close all;

%% Timing and visualization
DT = 0.2;          % control period [s]
ARROW_SCALE = 0.4; % quiver length for headings (meters)

%% Vicon setup (reuse template style)
VICON_CLIENT = Vicon.Client();
VICON_CLIENT.destroy();
VICON_CLIENT.initialize();

ERLICH_VICON  = "Object1";
ERLICH_PORT   = 50804;
BACHMAN_VICON = "Object2";
BACHMAN_PORT  = 50805;

bot1 = turtlebot(VICON_CLIENT, ERLICH_VICON, ERLICH_PORT);
bot2 = turtlebot(VICON_CLIENT, BACHMAN_VICON,  BACHMAN_PORT);

%% Controller parameters (meters)
% Per Example 5.2 in Sanfelice et al. (ACC 2006):
% - Single target at (3, 0)
% - Obstacle at (1, 0) with radius 1/(20*sqrt(2))
% - Mode q determines PATH around obstacle, not destination

% Single target (both robots go here via different paths)
params.xt = 3.0;
params.yt = 0.0;

% Obstacle parameters (matching paper)
params.obs_cx = 1.0;
params.obs_cy = 0.0;
params.obs_R  = 1/(20*sqrt(2));  % δ = 1/(20*sqrt(2)) ≈ 0.0354 from paper

% Wedge angle for mode regions (per Figure 5)
params.wedge_angle = pi/4;  % 45 degrees

hyb.mu     = 1.2;
hyb.lambda = 0.05;
hyb.gamma  = 0.01;

% Mapping gains from Cartesian velocity to unicycle commands
k_v     = 0.5;
k_theta = 1.5;

noise_std = 0.0; % measurement noise is already present via Vicon

% Initial modes select distinct classes around the obstacle
q_bot1 = 1; % takes upper path around obstacle
q_bot2 = 2; % takes lower path around obstacle

t = 0;
recorded_data = [];

%% Main control loop
while true
    % 1) Read odometry (millimeters from Vicon)
    [m_x_mm, m_y_mm, m_theta, bot1] = bot1.odom();
    [c_x_mm, c_y_mm, c_theta, bot2] = bot2.odom();

    % 2) Safety bounds copied from client_circle template
    if m_y_mm > 2400 || m_y_mm < -2700 || m_x_mm > 1900 || m_x_mm < -2000 || ...
       c_y_mm > 2400 || c_y_mm < -2700 || c_x_mm > 1900 || c_x_mm < -2000
        bot1.stop_car();
        bot2.stop_car();
        disp("out of bounds");
        break;
    end

    % 3) Convert to meters for the hybrid controller
    m_state = [m_x_mm; m_y_mm] / 1000;
    c_state = [c_x_mm; c_y_mm] / 1000;

    % 4) Hybrid controller using measurement y for both switching and control
    [~, q_bot1_next, info_m] = hybrid_controller_step(m_state, q_bot1, params, hyb, DT, noise_std);
    [~, q_bot2_next, info_c] = hybrid_controller_step(c_state, q_bot2, params, hyb, DT, noise_std);

    % Desired planar velocities (already based on measurement)
    u_m = info_m.u;
    u_c = info_c.u;

    % 5) Map to unicycle commands (heading alignment)
    psi_m = atan2(u_m(2), u_m(1));
    heading_error_m = wrapToPi(psi_m - m_theta);
    m_v_cmd     = min(max(k_v * norm(u_m), bot1.min_v), bot1.max_v);
    m_gamma_cmd = min(max(k_theta * heading_error_m, bot1.min_gamma), bot1.max_gamma);

    psi_c = atan2(u_c(2), u_c(1));
    heading_error_c = wrapToPi(psi_c - c_theta);
    c_v_cmd     = min(max(k_v * norm(u_c), bot2.min_v), bot2.max_v);
    c_gamma_cmd = min(max(k_theta * heading_error_c, bot2.min_gamma), bot2.max_gamma);

    % 6) Send commands
    bot1 = bot1.drive(m_v_cmd, m_gamma_cmd, DT);
    bot2 = bot2.drive(c_v_cmd, c_gamma_cmd, DT);

    q_bot1 = q_bot1_next;
    q_bot2 = q_bot2_next;

    % 7) Log data
    recorded_data = [recorded_data; ...
        t, m_state.', m_theta, q_bot1, c_state.', c_theta, q_bot2, ...
        m_v_cmd, m_gamma_cmd, c_v_cmd, c_gamma_cmd];

    % 8) Visualization for quick monitoring
    figure(1); clf; hold on; axis equal;
    th = linspace(0, 2*pi, 200);
    plot(params.obs_cx + params.obs_R*cos(th), params.obs_cy + params.obs_R*sin(th), 'k-', 'LineWidth', 2);
    plot(params.xt, params.yt, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    plot(m_state(1), m_state(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    quiver(m_state(1), m_state(2), ARROW_SCALE*cos(m_theta), ARROW_SCALE*sin(m_theta), 'b');
    plot(c_state(1), c_state(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    quiver(c_state(1), c_state(2), ARROW_SCALE*cos(c_theta), ARROW_SCALE*sin(c_theta), 'g');
    title(sprintf('t = %.1f s', t));
    xlabel('x [m]'); ylabel('y [m]');
    drawnow;

    % 9) Goal check (both robots reach epsilon ball)
    eps_goal = 0.1;
    goal = [params.xt; params.yt];
    if norm(m_state - goal) < eps_goal && ...
       norm(c_state - goal) < eps_goal
        disp('Both robots reached the goal.');
        break;
    end

    pause(DT/2);
    t = t + DT;
end

% Stop both robots at exit
bot1.drive(0, 0, 0);
bot2.drive(0, 0, 0);

% Post-processing plots (optional refinement)
if ~isempty(recorded_data)
    figure; plot(recorded_data(:,1), recorded_data(:,5), 'b-', 'LineWidth', 1.5); hold on;
    plot(recorded_data(:,1), recorded_data(:,8), 'g-', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Mode q'); legend('Bot1', 'Bot2'); grid on;
    title('Mode evolution during experiment');

    goal = [params.xt; params.yt];
    dist_bot1 = vecnorm(recorded_data(:,2:3).' - goal);
    dist_bot2 = vecnorm(recorded_data(:,6:7).' - goal);
    figure; plot(recorded_data(:,1), dist_bot1, 'b-', 'LineWidth', 1.5); hold on;
    plot(recorded_data(:,1), dist_bot2, 'g-', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Distance to goal [m]'); legend('Bot1', 'Bot2'); grid on;
    title('Convergence to goal');
end
