function [u_opt] = MPC_controller(input)
    arguments
        input (4,1) double {mustBeReal}  
    end
    
    x_current = input;  
    
    % System matrices 
    Ad = [1.002 0 0 0; 
          0 1.002 0 0; 
          0 0.1103 0.9909 0.0004; 
          0 0.3372 -0.0090 0.9989];
    Bd = [0; 0; 0.0412; 0.0407];
    
    % MPC parameters
    N = 10;  % Reduced for speed
    Q = diag([100, 1, 0.1, 0.1]);  % 4x4 state weights
    R = 1;   % input weight
    P = 10 * Q;
    u_max = 5;
    du_max = 2;
    x_ref = zeros(4,1);
    
    % Setup optimization
    n = size(Ad, 1);  % 4 states
    m = size(Bd, 2);  % 1 input
    
    u = sdpvar(m, N, 'full');
    x = sdpvar(n, N+1, 'full');
    
    % Constraints
    constraints = [x(:,1) == x_current];
    
    for k = 1:N
        constraints = [constraints, x(:,k+1) == Ad*x(:,k) + Bd*u(:,k)];
    end
    
    % Input magnitude constraints
    constraints = [constraints, -u_max <= u <= u_max];
    
    % Input rate constraints (simplified)
    for k = 2:N
        du = u(:,k) - u(:,k-1);
        constraints = [constraints, -du_max <= du <= du_max];
    end
    
    % Objective
    objective = 0;
    for k = 1:N
        objective = objective + (x(:,k) - x_ref)' * Q * (x(:,k) - x_ref);
        objective = objective + u(:,k)' * R * u(:,k);
    end
    objective = objective + (x(:,N+1) - x_ref)' * P * (x(:,N+1) - x_ref);
    
    % Solve
    diagnostics = optimize(constraints, objective);
    
    if diagnostics.problem == 0
        u_opt = value(u(:,1));
    else
        u_opt = 0;
        warning('MPC infeasible');
    end
end