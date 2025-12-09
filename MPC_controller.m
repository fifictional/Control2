function [u_opt] = MPC_controller(input)
    arguments
        input (5,:) = [0,0]
    end
    x_current = input(1);
    u_prev = input(2);


    Ad = [1.002 0 0 0; 0 1.002 0 0; 0 0.1103 0.9909 0.0004; 0 0.3372 -0.0090 0.9989];
    Bd = [0;
         0;
    0.0412;
    0.0407];
    
    N = 50; % prediction horizon 
    Q = [100   0     0     0;  
       0    1     0     0;  
       0    0   0.1   0;     
       0    0     0   0.1];  
    % state weightinig matrix
    R = 1;   % input weighting matrix                   
    P = 10 * Q;    % terminal cost  
    u_max = 5;
    du_max = 2;
    x_ref = [0; 0; 0; 0];

    n = size(Ad, 1);  % number of states
    m = size(Bd, 2);  % number of inputs
    
    u = sdpvar(m, N, 'full');    % control sequence,  a symbolic optimisation variable in YALMIP
    x = sdpvar(n, N+1, 'full');  % state sequence
    disp(u)
    disp(x)
    
    constraints = [];
    objective = 0;
    
    constraints = [constraints, x(:,1) == x_current];
    
    for k = 1:N
        constraints = [constraints, x(:,k+1) == Ad*x(:,k) + Bd*u(:,k)];
    end
    
    constraints = [constraints, -u_max <= u <= u_max]; % input mag limit
    
    for k = 1:N % input rate limit
        if k == 1
            du = u(:,k) - u_prev;
        else
            du = u(:,k) - u(:,k-1);
        end
        constraints = [constraints, -du_max <= du <= du_max];
    end
    
    % objective function 
    for k = 1:N
        disp(x(:,k))
        % stage cost
        objective = objective + (x(:,k) - x_ref)' * Q * (x(:,k) - x_ref);
        objective = objective + u(:,k)' * R * u(:,k);
    end
    
    % terminal cost
    objective = objective + (x(:,N+1) - x_ref)' * P * (x(:,N+1) - x_ref);

    
    diagnostics = optimize(constraints, objective);
    
    if diagnostics.problem == 0
        u_opt = value(u(:,1));      
        feasible = true;
        cost = value(objective);
    else
        u_opt = 0;                  
        feasible = false;
        cost = inf;
        warning('MPC infeasible');
    end
end
