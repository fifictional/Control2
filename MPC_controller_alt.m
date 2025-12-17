function [u_opt] = MPC_controller(input)
%#codegen
    persistent optimizerObj u_prev
    
    x_current = input;  % [theta; alpha; theta_dot; alpha_dot]
    
    % Initialize persistent variables
    if isempty(optimizerObj)
        % Define system
        Ad = [1 0 0.002 0; 
              0 1 0 0.002; 
              0 0.1103 0.9909 0.0004; 
              0 0.3372 -0.0090 0.9989];

        Bd = [0; 0; 0.0412; 0.0407];
        
        N = 100;  % Reduced horizon for speed
        Q = diag([1, 1, 1, 2]); % 0.5, 0.5, 0.2, 0.3
        R = 0.0001; %0.001
        u_max = 10;
        du_max = 10;
        
        % Define optimization variables
        x0 = sdpvar(4,1);      % Initial state
        U = sdpvar(1,N);       % Control sequence
        
        % Build constraints and objective
        constraints = [];
        objective = 0;
        x = x0;
        
        for k = 1:N
            % Dynamics
            x = Ad*x + Bd*U(k);
            
            % Input magnitude constraints
            constraints = [constraints, -u_max <= U(k) <= u_max];
            
            % Input rate constraints (for k>1)
            if k > 1
                du = U(k) - u_prev;
                constraints = [constraints, -du_max <= du <= du_max];
            end
            
            % Cost
            objective = objective + x'*Q*x + U(k)'*R*U(k);
        end
        
        % Terminal cost
        [~, P, ~] = dlqr(Ad, Bd, Q, R);
        % Terminal cost
        objective = objective + x'*P*x;
        
        % Create OPTIMIZER object (compiles once!)
        ops = sdpsettings('verbose', 0, 'solver', 'quadprog');
        optimizerObj = optimizer(constraints, objective, ops, x0, U(1));
        
        % Initialize previous control
        u_prev = 0;
    end
    
    % Solve MPC using pre-compiled optimizer (FAST!)
    [u_opt, errorcode] = optimizerObj(x_current);
    
    % Handle errors
    if errorcode ~= 0
        u_opt = 0;
        warning('MPC solver error');
    end
    
    % Update previous control
    u_prev = u_opt;
end