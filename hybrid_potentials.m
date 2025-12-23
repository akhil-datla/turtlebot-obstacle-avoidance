function [V1, V2, gradV1, gradV2] = hybrid_potentials(x, params)
% hybrid_potentials  Mode-dependent Lyapunov functions/gradients.
%   Implements the potentials from Sanfelice et al. (ACC 2006), Eq. (11):
%       V_q(x) = 0.5*||x - a||^2 + B(d_q(x))
%   where d_q(x) is the distance to the complement of the free set O_q and
%   B(·) is the logarithmic barrier keeping trajectories away from the obstacle
%   and switching manifold.
%
%   Per Figure 5 in the paper, the two-line wedge geometry:
%       - Two lines at ±θ through obstacle center form an X
%       - Upper wedge: region above BOTH lines
%       - Lower wedge: region below BOTH lines
%
%   Mode definitions (matching paper Figure 5):
%       q = 1 : O_1 = excludes upper wedge → robot takes LOWER path
%       q = 2 : O_2 = excludes lower wedge → robot takes UPPER path
%
%   Per Example 5.2, there is a single target.
%
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
%   Implements the TWO-LINE wedge geometry from Sanfelice et al. (ACC 2006)
%   Figure 5. Two lines at ±θ through the obstacle center form an X:
%       Line A: y = cy + tan(θ)*(x - cx)  [positive slope]
%       Line B: y = cy - tan(θ)*(x - cx)  [negative slope]
%
%   Wedge regions:
%       Upper wedge: points ABOVE both lines (y > both line values)
%       Lower wedge: points BELOW both lines (y < both line values)
%
%   Mode definitions (per paper Figure 5):
%       Mode 1 (O_1): Excludes UPPER wedge → robot takes LOWER path
%       Mode 2 (O_2): Excludes LOWER wedge → robot takes UPPER path
%
%   d_q(x) is the distance to the COMPLEMENT of O_q (the excluded wedge):
%       - If x is in O_q (valid): d_q = distance to wedge boundary (positive)
%       - If x is in complement (forbidden wedge): d_q ≈ 0, barrier → ∞

% Minimum distance threshold to avoid singularities in barrier function
MIN_DIST = 1e-6;

% Get wedge angle (default 45 degrees per paper)
if isfield(params, 'wedge_angle')
    theta_wedge = params.wedge_angle;
else
    theta_wedge = pi/4;
end

cx = params.obs_cx;
cy = params.obs_cy;
R  = params.obs_R;

% Position relative to obstacle center
dx = x(1) - cx;
dy = x(2) - cy;

% Distance to circular obstacle boundary
r = sqrt(dx^2 + dy^2);
if r < MIN_DIST
    r = MIN_DIST;
end
dist_obs = r - R;
grad_obs = [dx; dy] / r;

% For points past the obstacle (x1 > cx + R), only obstacle barrier matters
if x(1) > cx + R
    d = max(dist_obs, MIN_DIST);
    grad_d = grad_obs;
    return;
end

% Two-line wedge geometry
slope = tan(theta_wedge);
norm_factor = sqrt(slope^2 + 1);

% Signed distances to each line:
%   Line A (positive slope): slope*dx - dy = 0
%   Line B (negative slope): -slope*dx - dy = 0
%
% Sign convention:
%   signed_dist_A > 0 means BELOW Line A
%   signed_dist_A < 0 means ABOVE Line A
%   signed_dist_B > 0 means BELOW Line B
%   signed_dist_B < 0 means ABOVE Line B
signed_dist_A = (slope * dx - dy) / norm_factor;
signed_dist_B = (-slope * dx - dy) / norm_factor;

switch q
    case 1
        % MODE 1: Exclude UPPER wedge → robot takes LOWER path
        % Upper wedge = ABOVE both lines = signed_dist_A < 0 AND signed_dist_B < 0

        in_upper_wedge = (signed_dist_A < 0) && (signed_dist_B < 0);

        if in_upper_wedge
            % Inside forbidden upper wedge - barrier should be very high
            dist_wedge = MIN_DIST;
            % Gradient points OUT of wedge (downward toward valid region)
            % Choose the line we're closer to
            if abs(signed_dist_A) < abs(signed_dist_B)
                grad_wedge = [slope; -1] / norm_factor;   % Normal to Line A, pointing down
            else
                grad_wedge = [-slope; -1] / norm_factor;  % Normal to Line B, pointing down
            end
        else
            % Outside upper wedge (in valid region O_1)
            % Distance to wedge = distance to the relevant boundary line
            % The boundary depends on which side of the obstacle center we're on
            if dx >= 0
                % Right of center: Line A forms the upper wedge boundary
                dist_wedge = signed_dist_A;  % Positive when below Line A
                grad_wedge = [slope; -1] / norm_factor;
            else
                % Left of center: Line B forms the upper wedge boundary
                dist_wedge = signed_dist_B;  % Positive when below Line B
                grad_wedge = [-slope; -1] / norm_factor;
            end
            dist_wedge = max(dist_wedge, MIN_DIST);
        end

    case 2
        % MODE 2: Exclude LOWER wedge → robot takes UPPER path
        % Lower wedge = BELOW both lines = signed_dist_A > 0 AND signed_dist_B > 0

        in_lower_wedge = (signed_dist_A > 0) && (signed_dist_B > 0);

        if in_lower_wedge
            % Inside forbidden lower wedge - barrier should be very high
            dist_wedge = MIN_DIST;
            % Gradient points OUT of wedge (upward toward valid region)
            if signed_dist_A < signed_dist_B
                grad_wedge = [-slope; 1] / norm_factor;   % Normal to Line A, pointing up
            else
                grad_wedge = [slope; 1] / norm_factor;    % Normal to Line B, pointing up
            end
        else
            % Outside lower wedge (in valid region O_2)
            % Distance to wedge = distance to the relevant boundary line
            if dx >= 0
                % Right of center: Line A forms the lower wedge boundary
                dist_wedge = -signed_dist_A;  % Positive when above Line A
                grad_wedge = [-slope; 1] / norm_factor;
            else
                % Left of center: Line B forms the lower wedge boundary
                dist_wedge = -signed_dist_B;  % Positive when above Line B
                grad_wedge = [slope; 1] / norm_factor;
            end
            dist_wedge = max(dist_wedge, MIN_DIST);
        end

    otherwise
        error('Unknown mode q=%d', q);
end

% Combine obstacle and wedge distances - use minimum (closest constraint)
if dist_obs <= 0
    % Inside obstacle - maximum penalty
    d = MIN_DIST;
    grad_d = grad_obs;
elseif dist_obs <= dist_wedge
    % Obstacle boundary is closer
    d = max(dist_obs, MIN_DIST);
    grad_d = grad_obs;
else
    % Wedge boundary is closer
    d = dist_wedge;
    grad_d = grad_wedge;
end
end
