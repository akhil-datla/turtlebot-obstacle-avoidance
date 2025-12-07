function [x_next, q_next, info] = hybrid_controller_step(x, q, params, hyb, dt, noise_std)
% hybrid_controller_step  One hybrid control step for the point-mass model.
%   Implements the flow/jump logic from Sanfelice et al. (ACC 2006),
%   for a two-mode controller. Mode q determines the path
%   around the obstacle (above or below), not the destination.
%   All controller computations use the measured state y = x + e.
%
% Inputs:
%   x         - current true state [2x1]
%   q         - current discrete mode (1 or 2)
%   params    - geometry and potential parameters
%   hyb       - struct with fields mu (>1), lambda (>0), gamma (>0)
%   dt        - integration step [s]
%   noise_std - measurement noise standard deviation [m]
%
% Outputs:
%   x_next  - propagated continuous state after one Euler step
%   q_next  - updated discrete mode via hysteresis logic
%   info    - struct logging potentials, control, measurement, and jump events

% 1) Measurement with optional Gaussian noise
e = zeros(2,1);
if noise_std > 0
    e = noise_std * randn(2,1);
end
y = x + e;

% 2) Potentials and gradients at the measurement (single target)
[V1_y, V2_y, gradV1_y, gradV2_y] = hybrid_potentials(y, params);

V_active = [V1_y; V2_y];
Vmin_y = min(V_active);
Vq_y = V_active(q);

% 3) Flow and jump conditions (Eqs. 12-16)
can_flow = (Vq_y <= hyb.mu * Vmin_y) || (Vq_y <= hyb.gamma);
can_jump = (Vq_y >= (hyb.mu - hyb.lambda) * Vmin_y);

% Prioritize jumps when both sets overlap
jumped = false;
q_new = q;
if can_jump
    [~, q_candidate] = min(V_active);
    if q_candidate ~= q
        q_new = q_candidate;
        jumped = true;
    end
end

% 4) Continuous control using the gradient at the measurement
if q_new == 1
    grad_active = gradV1_y;
else
    grad_active = gradV2_y;
end
u = -grad_active;

% 5) Euler integration of the true continuous state
x_next = x + dt * u;
q_next = q_new;

% 6) Log auxiliary data
info = struct();
info.V_y = V_active;
info.q_before = q;
info.q_after = q_next;
info.jumped = jumped;
info.u = u;
info.y = y;
info.e = e;
end
