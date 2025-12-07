function dBdz = barrier_derivative_numeric(d)
% barrier_derivative_numeric  Analytical derivative of the barrier potential.
%   Despite the historical name, this function implements the exact derivative
%   dB/dz from Sanfelice et al. (ACC 2006) for
%       B(z) = (z-1)^2 * log(1/z), 0 < z < 1.
%   The derivative is
%       dB/dz = -2*(z-1)*log(z) - (z-1)^2 / z.
%   For z <= 0 the derivative is undefined and is returned as 0 to avoid
%   corrupting the gradient when the state is inside the obstacle. For z >= 1
%   the barrier is flat and the derivative is 0.

if d <= 0 || d >= 1
    dBdz = 0;
    return;
end

% Clamp to avoid log(0) without changing the analytical structure.
eps_d = 1e-6;
z = min(max(d, eps_d), 1 - eps_d);

% Exact derivative of the barrier.
dBdz = -2 * (z - 1) .* log(z) - ((z - 1).^2) ./ z;
end
