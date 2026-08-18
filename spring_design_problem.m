function [f, c, ceq] = spring_design_problem(x)
%SPRING_DESIGN_PROBLEM  Standard tension/compression coil spring design problem (P2).
%
%   Minimize the total weight of a coil spring.
%
%   x = [x1 x2 x3] = [d D P]
%        d: wire diameter, D: mean coil diameter, P: number of active coils
%
%   Reference formulation:
%     Belegundu A D, Arora J S. A study of mathematical programming methods
%       for structural optimization. Part I/II. IJNME, 1985, 21(9).
%     Cagnina L C, Esquivel S C, Coello C A C. Solving engineering
%       optimization problems with the simple constrained particle swarm
%       optimizer. Informatica, 2008, 32(3):319-326.
%     (classic optimum f* = 0.012665 at
%      x* = (0.051690, 0.356750, 11.287126))
%
%   Outputs:
%     f    objective value (weight)
%     c    inequality constraints, c(x) <= 0
%     ceq  equality constraints (none)

    x1 = x(1);  x2 = x(2);  x3 = x(3);

    % ------------------------------------------------------------------
    % Objective
    % ------------------------------------------------------------------
    f = x1^2*x2*(x3 + 2);

    % ------------------------------------------------------------------
    % Inequality constraints  c(x) <= 0   (standard 4-constraint form)
    % ------------------------------------------------------------------
    c = [ 1 - x2^3*x3/(71785*x1^4);                                    % g1: shear stress
          (4*x2^2 - x1*x2)/(12566*(x2*x1^3 - x1^4)) + 1/(5108*x1^2) - 1; % g2: shear/wahl factor
          1 - 140.45*x1/(x2^2*x3);                                     % g3: surge frequency
          (x1 + x2)/1.5 - 1 ];                                         % g4: deflection
    ceq = [];
end
