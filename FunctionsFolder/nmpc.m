function [x, u] = nmpc(runningcosts, ...
              constraints, ...
              linearconstraints, system, ...
              mpciterations, N, T,  xmeasure, u0, x_ref, x_TV_measure, x_TV_ref, ...
              Determine_X_new, EV_system,TV_system, computeCurvedBicycleModel, ...
              printHeader, printClosedloopData, plotTrajectories ,TV_prediction, TV_dynamics)

%the functions ARE:
%   runningcosts:         evaluates the running costs for state and control
%                         at one sampling instant.
%                         The function returns the running costs for one
%                         sampling instant.
%          Usage: [cost] = runningcosts(x, u, x_ref,EV_system)
%                 with  state x, reference state x_ref, control u and
%                 system EV_systemeters EV_system
%   constraints:        
%                         The function returns the value of the
%                         restrictions for a sampling instance separated
%                         for inequality restrictions c and equality
%                         restrictions ceq.
%          
%   linearconstraints:    sets the linear constraints of the discretized
%                         optimal control problem. This is particularly
%                         useful to set control and state bounds.
%                         The function returns the required matrices for
%                         the linear inequality and equality constraints A
%                         and Aeq, the corresponding right hand sides b and
%                         beq as well as the lower and upper bound of the
%                         control.
%          Usage: [A, b, Aeq, beq, lb, ub] = linearconstraints()
%
%   system:               evaluates the linearized dynamics equation describing the
%                         process given time scale T, state vector x and control
%                         u, the linearisation state x_curr and the
%                         dynamics matrices linearisation
%                         The function returns the state vector x at the
%                         next time instant.
%          
%   Determine_X_new:      evaluates the dynamics equation
%                         given  state vector x and control  u.
%                         The function returns the state vector x at the
%                         next time instant.
%          
%   computeCurvedBicycleModel:  evaluates the matrices for the linearized dynamics
%                         given  state vector x and control  u and time scale T.
%                         
%          Usage: [A_d, B_d, x0cft] = computeCurvedBicycleModel(x, u, T, EV_system)
%                 with time scale T, state x and control u, EV_system is the system
%                 Properties
%   plotTrajectories:    Graphical output of the trajectories (Visualisation)
%
% Arguments are:
%   mpciterations:  Number of MPC iterations to be performed
%   N:              Length of optimization horizon
%   T:              Sampling interval
%   xmeasure:       State measurement of initial value of EV
%   u0:             Initial guess of open loop control
%   x_ref:          Reference trajectory for the EV
%   EV_system: A structure  that stores the propterties of the Problem

    tol_opt = 1e-6;
    % options for fmincon
    options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 10000,...
                'Algorithm', 'active-set',...
                'FinDiffType', 'forward',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1,...
                'TolCon', 1e-6, ...
                'TolConSQP', 1e-6); 
            

    % Initilize closed loop data
    x = [];  % Position of the EV
    u = [];
    cost_array = 0;
    cost = 0;   % the running Cost
    Cost_total = cost; % total Cost of the whole Simulation
    
    
    x_TV_predicted = []; % position of tv car
    w = TV_system.noise;
    
    % input of the previous time step
    u_last = [0;0];
    
    % the NMPC iteration FROM THE Book
    mpciter = 0;
    while(mpciter < mpciterations)
        % Step 1:
        % Obtain new initial value
        x0 = xmeasure;       
        
        % Obtain matrices from Curved Bicycle Model
        [Ad, Bd, x0fcT] = computeCurvedBicycleModel(x0, T, 0, EV_system.l_f, EV_system.l_r);
        myMatrices.A=Ad;
        myMatrices.B= Bd;
        myMatrices.x0fcT = x0fcT;
        
        % Predict Trajectory of TV  for N horizon
        x_TV_0 = x_TV_measure; 
        x_TV_ref_N = x_TV_ref(:,mpciter+1:mpciter+N+1);
        
        X_TV = TV_prediction(x_TV_0,x_TV_ref_N,N,T,w);
        
        % Update TV position  : accumulate real Status of TV from the
        % Beginning
        x_TV_predicted = [ x_TV_predicted, x_TV_measure ];
        
        % TV next status (with noise)
        [x_TV_measure,~] = TV_dynamics(x_TV_0, x_TV_ref(:,mpciter+1),T, w);
 
        
        % Step 2 :
        %   Solve the optimal control problem
        t_Start = tic;
        [u_new, V_current, exitflag, output] = solveOptimalControlProblem ...
            (runningcosts, constraints, ...
             linearconstraints, system, ...
            N, T, x0, u0, options, xmeasure, x_ref(:,mpciter+1:mpciter+N+1),X_TV, ...
        EV_system , TV_system, myMatrices, u_last,Determine_X_new);

        t_Elapsed = toc( t_Start );
        
        % the predicted trajectory of the EV
        x_pred = dynamics(system,N,T,x0,u_new, EV_system, myMatrices,Determine_X_new);
        
        %  Print solution
        printSolution(printHeader, printClosedloopData, ...
                      plotTrajectories, mpciter, x0, u_new, ...
                      exitflag, output, t_Elapsed, EV_system, x_ref, x_pred, x_TV_ref, x_TV_measure);


        %   Store closed loop data
        x = [ x xmeasure ]; 
        u = [ u u_new(:,1)];
        
        % Step 3:
        %  Apply control to process
        xmeasure = applyControl(system, T, x0, u_new, myMatrices);
        
        % add previous Cost value
        run_cost =runningcosts(xmeasure, [u_last,u_new(:,1)], x_ref(:,1), EV_system);
        cost_array = [cost_array run_cost];
        cost = cost + run_cost;
        Cost_total = [Cost_total cost];
 
        
% %         % Cost Plot
%           plotCost(Cost_total,cost_array, mpciter);
% % % 
% % %         % Velocity Graph
 %          plotVelocity(x, EV_system, mpciterations, mpciter);
% % % 
% % %         % Position d
  %         plot_d(mpciter, mpciterations,N, x);

        % plot Cars allure
        plotTrajectories( x, x0, x_ref,EV_system,x_TV_measure)
        
        % refresh input for next iter
        u0 = shiftHorizon(u_new);
        u_last = u_new(:,1);
   
        mpciter = mpciter+1  ;
    end
        
end

function xapplied = applyControl(system, T, x0, u, myMatrices)
    xapplied = system(x0, u(:,1), T, x0, myMatrices);
end


function [u, V, exitflag, output] = solveOptimalControlProblem ...
    (runningcosts, constraints, ...
    linearconstraints, system, N,T, x0, u0, ...
    options, xmeasure, x_ref,x_TV, EV_system,TV_system, myMatrices, u_last,Determine_X_new)
    % Set control and linear bounds
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    % update velocity depedent  part in a_r_k (TV  rectangle a)
    a_r = [];
    
    for k=1:N
        [Anew, bnew, Aeqnew, beqnew, lbnew, ubnew] = ...
               linearconstraints();
        A = blkdiag(A,Anew);
        b = [b, bnew];
        Aeq = blkdiag(Aeq,Aeqnew);
        beq = [beq, beqnew];
        lb = [lb, lbnew];
        ub = [ub, ubnew];
        a_min_longitudinal = -9 ; 
    
        a_r_k = -1/a_min_longitudinal * max([0,(xmeasure(4)^2- x_TV(2,k)^2)]);
        a_r = [a_r a_r_k];
    end

    % update a_r 
    TV_system.a_r_velocity_dependent = TV_system.a_r_velocity_independent + a_r;
    
    % Solve optimization problem
    [u, V, exitflag, output] = fmincon(@(u) costfunction(runningcosts, ...
    system, N,T, x0, ...
    u,x_ref,EV_system, myMatrices, u_last,Determine_X_new), u0, A, b, Aeq, beq, lb, ...
    ub, @(u) nonlinearconstraints(constraints, ...
    system, N,T, x0, u, EV_system, TV_system, x_ref,x_TV,myMatrices,Determine_X_new), options);
        

end

function x = computeOpenloopSolution(system, N, T, x0, u, myMatrices)
    x(:,1) = x0;
    for k=1:N
        x(:,k+1) = system(x(:,k), u(:,k), T, x0, myMatrices);                       
    end
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))];
end

function cost = costfunction(runningcosts, system, ...
                    N,T, x0, u, x_ref, EV_system, ...
                    myMatrices, u_last,Determine_X_new)
    cost = 0;
    n = length(x0);
    x = zeros(n, N+1);
    % x: openloopsolution
    x = dynamics(system,N,T,x0,u, EV_system, myMatrices,Determine_X_new); 
    cost = cost+ runningcosts(x(:,2), [u_last,u(:,1)], x_ref(:,2), EV_system);
    
    for k=2:N
        cost = cost+runningcosts(x(:,k+1), u(:,k-1:k), x_ref(:,k+1), EV_system);
    end
end

function [c,ceq] = nonlinearconstraints(constraints, ...
    system, ...
    N,T, x0, u, EV_system, TV_system, x_ref,x_TV,myMatrices, Determine_X_new)
    x = zeros(N+1, length(x0));

    %OpenloopSolution
    x = dynamics(system,N,T,x0,u, EV_system, myMatrices, Determine_X_new);
    c = [];
    ceq = [];
    
    for k=1:N
        [cnew, ceqnew] = constraints(x(:,k),x_ref(:,k),x_TV(:,k), k, EV_system, TV_system);
        c = [c cnew];
        ceq = [ceq ceqnew];
    end
    % terminal constraints
%     [cnew, ceqnew] = constraints(x(:,N+1),x_ref(:,N+1),x_TV(:,N+1), N+1, EV_system, TV_system);
%     c = [c cnew];
%     ceq = [ceq ceqnew];
end


%%%%%% OUTPUT %%%%%%%%%%%%%
function printSolution(printHeader, printClosedloopData, ...
                       plotTrajectories, mpciter, x0, u,...
                       exitflag, output, t_Elapsed, EV_system, x_ref, x , x_refTV, x_TV)
    % Print results
    if (mpciter == 0)
        printHeader();
    end
    
    printClosedloopData(mpciter, u, x0, t_Elapsed, x_TV);
    switch exitflag
        case -2
            fprintf(' Error: exitflag = -2 \n')
        
        case -1
            fprintf(' Error:exitflag = -1 \n')
        
        case 0
            fprintf(' no error : exitflag = 0 \n')
    end
    %plotTrajectories(x, x0, x_ref,EV_system,x_refTV ) 
end


function X = dynamics(system, N,T, x0, u0, EV_system, myMatrices,Determine_X_new)
%% Calculate the state vector X 

    X = x0;

    % Calculate  state vector
    xi = computeOpenloopSolution(system, N, T, x0, u0(:,1+N*(1-1):N*1), myMatrices);
    % update 
    x0 = xi(:,N+1); 
    xi = xi(:,2:N+1); % delete x0

    X = [X,xi];
     

end

%% Utils
function plotVelocity(x, EV_system, mpciterations, mpciter)
    figure(3);
    plot(0:mpciter,x(4,:),'r');
    axis([0 mpciterations 0 EV_system.v_max]);
    title('Velocity  NMPC TEST ');
    xlabel('time (k)');
    ylabel('velocity v (m/s)');
    saveas(figure(3),[pwd '/plots/second_try/velocity.jpg']);

end

function plotCost(Cost_total,cost_array, mpciter)
   figure(4)
   plot(0:mpciter,Cost_total(1:mpciter+1),'r');
   title('Standard MPC');
   xlabel('time instant k');
   ylabel(' total Cost function value');
   saveas(figure(4),[pwd '/plots/second_try/totalCost.jpg']);

   
    figure(9)
   plot(0:mpciter,cost_array(1:mpciter+1),'r');
   title('Standard MPC');
   xlabel('time instant k');
   ylabel('Cost function value');
   saveas(figure(9),[pwd '/plots/second_try/cost.jpg']);

   
end

function plot_d(mpciter, mpciterations,N, x)
    figure(5);
    plot(0:mpciter,x(2,:),'r');
    axis([0 mpciterations+N+1 -5 5]);
    title('d  NMPC TEST ');
    xlabel('time (k)');
    ylabel('d (m)');
    
    
    saveas(figure(5),[pwd '/plots/second_try/d.jpg']);
end


