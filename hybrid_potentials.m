function [V1, V2, gradV1, gradV2] = hybrid_potentials(x, params)
% hybrid_potentials  Mode-dependent Lyapunov functions/gradients.
%   Implements the potentials from Sanfelice et al. (ACC 2006), Eq. (11):
%       V_q(x) = 0.5*||x - a||^2 + B(d_q(x))
%   where d_q(x) is the distance to the complement of the free set O_q and
%   B(·) is the logarithmic barrier keeping trajectories away from the obstacle
%   and switching manifold. Modes correspond to two homotopy classes:
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
%   WEDGES emanating from the obstacle, not half-planes:
%       O_1 = points that can reach goal by going ABOVE the obstacle
%       O_2 = points that can reach goal by going BELOW the obstacle
%   
%   Both regions include the goal and all points past the obstacle.
%   The wedge angle determines how "wide" each region is.

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
% Line passes through (cx, cy) with slope ±tan(theta_wedge)
switch q
    case 1
        % O_1: below the upper wedge line (y - cy < tan(θ) * (x - cx))
        % Wedge line: y = cy + tan(θ) * (x - cx)
        % Distance to line ax + by + c = 0 where:
        %   tan(θ)*(x - cx) - (y - cy) = 0
        %   tan(θ)*x - y - tan(θ)*cx + cy = 0
        slope = tan(theta_wedge);
        % Signed distance (positive when below line, negative when above)
        signed_dist = (slope * (x(1) - cx) - (x(2) - cy)) / sqrt(slope^2 + 1);
        
        if signed_dist < 0
            % Above the line (in wrong region for mode 1)
            dist_wedge = abs(signed_dist);
            % Gradient points toward the line (downward and to the right)
            % Normal to line ax + by + c = 0 is [a; b], here [slope; -1]
            % We want to point toward valid region (below line), so negate
            grad_wedge = -[slope; -1] / sqrt(slope^2 + 1);
        else
            % Below the line (in valid region)
            dist_wedge = signed_dist;
            % Gradient points away from line (downward and to the right)
            grad_wedge = [slope; -1] / sqrt(slope^2 + 1);
        end
        
    case 2
        % O_2: above the lower wedge line (y - cy > -tan(θ) * (x - cx))
        slope = -tan(theta_wedge);
        signed_dist = ((x(2) - cy) - slope * (x(1) - cx)) / sqrt(slope^2 + 1);
        
        if signed_dist < 0
            % Below the line (in wrong region for mode 2)
            dist_wedge = abs(signed_dist);
            % Gradient points toward the line (upward and to the right)
            % We want to point toward valid region (above line), so negate
            grad_wedge = -[-slope; 1] / sqrt(slope^2 + 1);
        else
            % Above the line (in valid region)
            dist_wedge = signed_dist;
            % Gradient points away from line (upward and to the right)
            grad_wedge = [-slope; 1] / sqrt(slope^2 + 1);
        end
        
    otherwise
        error('Unknown mode q=%d', q);
end

% Combine obstacle and wedge distances
if dist_obs <= 0
    % Inside obstacle
    d = MIN_DIST;
    grad_d = grad_obs;
elseif dist_obs <= dist_wedge
    d = max(dist_obs, MIN_DIST);
    grad_d = grad_obs;
else
    d = max(dist_wedge, MIN_DIST);
    grad_d = grad_wedge;
end
end
