function B = barrier_B(d)
% barrier_B  Barrier potential used in the hybrid Lyapunov functions.
%   B(d) implements the logarithmic barrier from Sanfelice et al. (ACC 2006),
%   Section IV, Eq. (11). The function is defined for the scaled distance
%   argument d \in \mathbb{R} as:
%       B(d) = (d-1)^2 * log(1/d),  for 0 < d < 1
%       B(d) = 0,                   for d >= 1
%       B(d) = +infinity,           for d <= 0  (implemented as a large constant)
%   Inputs are clipped to avoid numerical issues near the singularity at d=0.
%
%   This barrier penalizes proximity to the obstacle boundary while remaining
%   identically zero when the state is at least one distance unit away from the
%   boundary. In the hybrid controller, the distance argument d_q(x) is the
%   signed distance to the complement of the mode-dependent free set O_q.

% Large constant used to approximate +infinity inside the obstacle.
INF_PENALTY = 1e6;

if d <= 0
    B = INF_PENALTY;
elseif d < 1
    % Clamp to avoid log(0) while preserving the analytical shape.
    eps_d = 1e-6;
    z = min(max(d, eps_d), 1 - eps_d);
    B = (z - 1).^2 .* log(1 ./ z);
else
    B = 0;
end
end
