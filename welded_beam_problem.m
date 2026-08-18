function [f, c, ceq] = welded_beam_problem(x)
%WELDED_BEAM_PROBLEM  Standard welded beam design problem (P1).
%
%   Minimize the fabrication cost of a welded beam.
%
%   x = [x1 x2 x3 x4] = [h l t b]
%        h: weld height, l: weld length, t: beam height, b: beam width
%
%   Reference formulation (standard 6-constraint version):
%     "A cost-effective algorithm for the solution of engineering problems
%      with particle swarm optimization", Engineering Optimization, 2010,
%      42(5):417-440, DOI:10.1080/03052150903305476.
%     (origin: Deb 1991; classic optimum f* = 1.724852 at
%      x* = (0.20573, 3.470489, 9.036624, 0.20573))
%
%   Outputs:
%     f    objective value (fabrication cost)
%     c    inequality constraints, c(x) <= 0
%     ceq  equality constraints (none)

    % ------------------------------------------------------------------
    % Fixed parameters
    % ------------------------------------------------------------------
    P         = 6000;      % load (lb)
    L         = 14;        % beam length (in)
    E         = 30e6;      % Young's modulus (psi)
    G         = 12e6;      % shear modulus (psi)
    tau_max   = 13600;     % allowable shear stress  (psi)
    sigma_max = 30000;     % allowable bending stress (psi)
    delta_max = 0.25;      % allowable end deflection (in)

    x1 = x(1);  x2 = x(2);  x3 = x(3);  x4 = x(4);

    % ------------------------------------------------------------------
    % Objective
    % ------------------------------------------------------------------
    f = 1.10471*x1^2*x2 + 0.04811*x3*x4*(14.0 + x2);

    % ------------------------------------------------------------------
    % Shear stress  tau(x)
    % ------------------------------------------------------------------
    R    = sqrt(x2^2/4 + ((x1 + x3)/2)^2);       % distance to weld centroid
    M    = P*(L + x2/2);                         % bending moment on weld
    J    = 2*sqrt(2)*x1*x2*(x2^2/12 + (x1+x3)^2/4);  % polar moment of weld group
    tau1 = P/(sqrt(2)*x1*x2);                    % primary shear stress
    tau2 = M*R/J;                                % torsional shear stress
    tau  = sqrt(tau1^2 + 2*tau1*tau2*(x2/(2*R)) + tau2^2);

    % ------------------------------------------------------------------
    % Bending stress, end deflection, buckling load
    % ------------------------------------------------------------------
    sigma = 6*P*L/(x4*x3^2);
    delta = 4*P*L^3/(E*x3^3*x4);
    Pc    = 4.013*E*sqrt(x3^2*x4^6/36)/L^2 * (1 - x3/(2*L)*sqrt(E/(4*G)));

    % ------------------------------------------------------------------
    % Inequality constraints  c(x) <= 0   (standard 6-constraint form)
    % ------------------------------------------------------------------
    c = [ tau   - tau_max;                                  % g1: shear stress
          sigma - sigma_max;                                % g2: bending stress
          x1    - x4;                                       % g3: h <= b
          0.10471*x1^2 + 0.04811*x3*x4*(14.0 + x2) - 5.0;   % g4
          delta - delta_max;                                % g5: end deflection
          P     - Pc ];                                     % g6: buckling load
    ceq = [];
end
