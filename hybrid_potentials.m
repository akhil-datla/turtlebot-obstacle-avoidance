function [V1, V2, gradV1, gradV2] = hybrid_potentials(x, params)
% hybrid_potentials  Mode-dependent Lyapunov functions/gradients.
%   Implements the potentials from Sanfelice et al. (ACC 2006), Eq. (11):
%       V_q(x) = 0.5*||x - a||^2 + B(d_q(x))
%   where d_q(x) is the distance to the complement of the free set O_q and
%   B(·) is the logarithmic barrier keeping trajectories away from the obstacle
%   and switching manifold. Modes correspond to two classes:
%       q = 1 : O_1 = wedge region for going ABOVE the obstacle
%       q = 2 : O_2 = wedge region for going BELOW the obstacle
%   Per Figure 5 in the paper, the regions are wedges emanating from the
%   obstacle, not half-planes. Both regions include the goal.
%   Per Example 5.2, there is a single target.
%   Inputs:
%       x        - 2x1 position vector [m]
%       params   - struct with fields obs_cx, obs_cy, obs_R, wedge_angle,
%                  and target position (params.xt/params.yt preferred,
%                  params.targets supported for legacy compatibility)
%   Outputs:
%       V1, V2     - potentials for modes 1 and 2 evaluated at x
%       gradV1, gradV2 - corresponding gradients

% Resolve the active goal
p_target = select_goal(params);

% Mode 1 distance and gradient
[d1, grad_d1] = distance_to_mode_complement(x, 1, params);

% Mode 2 distance and gradient
[d2, grad_d2] = distance_to_mode_complement(x, 2, params);

% Base quadratic term
base_grad = x - p_target;
base_val = 0.5 * sum((x - p_target).^2);

% Barrier contributions
B1 = barrier_B(d1);
B2 = barrier_B(d2);

% Barrier gradients
B1_grad = barrier_derivative_numeric(d1) .* grad_d1;
B2_grad = barrier_derivative_numeric(d2) .* grad_d2;

% Potentials and gradients
V1 = base_val + B1;
V2 = base_val + B2;

gradV1 = base_grad + B1_grad;
gradV2 = base_grad + B2_grad;
end

function goal = select_goal(params)
% select_goal  Returns the 2x1 goal vector from params.
%   Per Example 5.2 in Sanfelice et al. (ACC 2006), there is a single target.

if isfield(params, 'xt') && isfield(params, 'yt')
    goal = [params.xt; params.yt];
elseif isfield(params, 'targets')
    % Legacy support: use first column of targets matrix
    goal = params.targets(:, 1);
else
    error('params must have xt/yt or targets field');
end
end

function [d, grad_d] = distance_to_mode_complement(x, q, params)
% distance_to_mode_complement  Computes d_q(x) and its gradient.
%
%   Per Sanfelice et al. (ACC 2006) Figure 5, the regions O_1 and O_2 are
%   defined by two lines through the obstacle center (cx, cy):
%       Line A: y = cy + tan(θ)*(x - cx)  [positive slope, at y=-1 when x=0]
%       Line B: y = cy - tan(θ)*(x - cx)  [negative slope, at y=+1 when x=0]
%
%   Mode regions (for x < cx, i.e., approaching the obstacle):
%       O_1 = points ABOVE Line A (upper path around obstacle)
%       O_2 = points BELOW Line B (lower path around obstacle)
%
%   At x=0 with cx=1, cy=0, θ=π/4:
%       O_1: y > -1  (includes robot at y=0.5)
%       O_2: y < +1  (includes robot at y=-0.5)
%       Overlap: -1 < y < 1 (synergistic region)
%
%   d_q(x) is the distance to the COMPLEMENT of O_q:
%       - If x is inside O_q: d_q = distance to boundary (positive)
%       - If x is outside O_q (in complement): d_q ≈ 0, triggering barrier

% Minimum distance threshold to avoid singularities in barrier function
MIN_DIST = 1e-6;

% Get wedge angle (default 45 degrees)
if isfield(params, 'wedge_angle')
    theta_wedge = params.wedge_angle;
else
    theta_wedge = pi/4; % 45 degrees
end

cx = params.obs_cx;
cy = params.obs_cy;
R  = params.obs_R;

% Distance to obstacle boundary (positive when outside, negative when inside)
r_vec = [x(1) - cx; x(2) - cy];
r = norm(r_vec);
dist_obs = r - R;

% Guard against degenerate radius
if r < 1e-9
    r = 1e-9;
    r_vec = [1e-9; 0]; % arbitrary direction
end

% Gradient of distance to obstacle (points away from center)
grad_obs = r_vec / r;

% For points past the obstacle (x1 > cx + R), both modes can reach goal
% Only apply wedge barrier for points that haven't cleared the obstacle
if x(1) > cx + R
    % Past the obstacle - only obstacle barrier matters
    d = max(dist_obs, MIN_DIST);
    grad_d = grad_obs;
    return;
end

% Compute distance to wedge boundary line
% Lines pass through (cx, cy) with slope ±tan(theta_wedge)
slope = tan(theta_wedge);
norm_factor = sqrt(slope^2 + 1);

switch q
    case 1
        % O_1: ABOVE Line A (for upper path around obstacle)
        % Line A: y = cy + tan(θ)*(x - cx)  [positive slope line]
        % At x=0 with cx=1: y = 0 + 1*(-1) = -1  (line is BELOW origin)
        % At x=2: y = 0 + 1*(1) = 1
        %
        % Line equation: tan(θ)*(x - cx) - (y - cy) = 0
        % Signed distance: POSITIVE when y < line (below, invalid for Mode 1)
        %                  NEGATIVE when y > line (above, valid for Mode 1)
        signed_dist = (slope * (x(1) - cx) - (x(2) - cy)) / norm_factor;

        if signed_dist < 0
            % Above Line A (in valid region O_1)
            dist_wedge = abs(signed_dist);
            % Gradient of d = |signed_dist| = -signed_dist
            % grad(signed_dist) = [slope; -1]/norm, so grad(d) = [-slope; 1]/norm
            % This points upward (away from line, into valid region)
            grad_wedge = [-slope; 1] / norm_factor;
        else
            % Below Line A (in complement of O_1 - WRONG region)
            dist_wedge = MIN_DIST;
            % Gradient points toward valid region (upward)
            grad_wedge = [-slope; 1] / norm_factor;
        end

    case 2
        % O_2: BELOW Line B (for lower path around obstacle)
        % Line B: y = cy - tan(θ)*(x - cx)  [negative slope line]
        % At x=0 with cx=1: y = 0 - 1*(-1) = 1  (line is ABOVE origin)
        % At x=2: y = 0 - 1*(1) = -1
        %
        % Line equation: tan(θ)*(x - cx) + (y - cy) = 0
        % Signed distance: POSITIVE when y > line (above, invalid for Mode 2)
        %                  NEGATIVE when y < line (below, valid for Mode 2)
        signed_dist = (slope * (x(1) - cx) + (x(2) - cy)) / norm_factor;

        if signed_dist < 0
            % Below Line B (in valid region O_2)
            dist_wedge = abs(signed_dist);
            % Gradient of d = |signed_dist| = -signed_dist
            % grad(signed_dist) = [slope; 1]/norm, so grad(d) = [-slope; -1]/norm
            % This points downward (away from line, into valid region)
            grad_wedge = [-slope; -1] / norm_factor;
        else
            % Above Line B (in complement of O_2 - WRONG region)
            dist_wedge = MIN_DIST;
            % Gradient points toward valid region (downward)
            grad_wedge = [-slope; -1] / norm_factor;
        end

    otherwise
        error('Unknown mode q=%d', q);
end

% Combine obstacle and wedge distances - take minimum (closest boundary)
if dist_obs <= 0
    % Inside obstacle - maximum penalty
    d = MIN_DIST;
    grad_d = grad_obs;
elseif dist_obs <= dist_wedge
    % Obstacle boundary is closer
    d = max(dist_obs, MIN_DIST);
    grad_d = grad_obs;
else
    % Wedge boundary is closer (or in wrong region where dist_wedge = MIN_DIST)
    d = max(dist_wedge, MIN_DIST);
    grad_d = grad_wedge;
end
end
