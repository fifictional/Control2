function [u_opt] = MPC_controller(input)
    %#codegen  
    
    arguments
        input (5,1) double {mustBeReal} 
    end
    
    % Extract states and previous control
    x_current = input(1:4);  
    u_prev = input(5);      
    
    % System matrices
    Ad = [1.002 0 0 0; 
          0 1.002 0 0; 
          0 0.1103 0.9909 0.0004; 
          0 0.3372 -0.0090 0.9989];
    Bd = [0; 0; 0.0412; 0.0407];
    
    % MPC parameters
    N = 50;
    Q = [100 0 0 0;  
         0 1 0 0;  
         0 0 0.1 0;     
         0 0 0 0.1];  
    R = 1;
    P = 10 * Q;
    u_max = 5;
    du_max = 2;
    x_ref = [0; 0; 0; 0];
    
    % Initialize output
    u_opt = 0; 
    
    % Setup optimization variables
    n = size(Ad, 1);
    m = size(Bd, 2);
    
    u = sdpvar(m, N, 'full');
    x = sdpvar(n, N+1, 'full');
    
    % Constraints
    constraints = [x(:,1) == x_current(:)];  
    
    for k = 1:N
        constraints = [constraints, x(:,k+1) == Ad*x(:,k) + Bd*u(:,k)];
    end
    
    constraints = [constraints, -u_max <= u <= u_max];
    
    % Rate constraints
    for k = 1:N
        if k == 1
            du = u(:,k) - u_prev;
        else
            du = u(:,k) - u(:,k-1);
        end
        constraints = [constraints, -du_max <= du <= du_max];
    end
    
    % Objective function
    objective = 0;
    for k = 1:N
        objective = objective + (x(:,k) - x_ref)' * Q * (x(:,k) - x_ref);
        objective = objective + u(:,k)' * R * u(:,k);
    end
    objective = objective + (x(:,N+1) - x_ref)' * P * (x(:,N+1) - x_ref);
    
    % Solve
    diagnostics = optimize(constraints, objective);
    
    % Output result
    if diagnostics.problem == 0
        u_opt = double(value(u(:,1)));
    else
        u_opt = double(0);
        warning('MPC infeasible');
    end
    
    % Ensure scalar output
    u_opt = u_opt(1);  
end