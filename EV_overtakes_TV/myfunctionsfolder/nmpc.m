function [x, u] = nmpc(runningcosts, ...
              constraints, ...
              linearconstraints, system, ...
              mpciterations, N, T, xmeasure, u0, x_ref, x_refTV, ...
              Determine_X_new, Properties_Obj, computeCurvedBicycleModel, ...
              printHeader, printClosedloopData, plotTrajectories)

%the functions ARE:
%   runningcosts:         evaluates the running costs for state and control
%                         at one sampling instant.
%                         The function returns the running costs for one
%                         sampling instant.
%          Usage: [cost] = runningcosts(x, u, x_ref,Properties_Obj)
%                 with  state x, reference state x_ref, control u and
%                 system Properties_Objeters Properties_Obj

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
%          Usage: [A_d, B_d, x0cft] = computeCurvedBicycleModel(x, u, T, Properties_Obj)
%                 with time scale T, state x and control u, Properties_Obj is the system
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
%   Properties_Obj: A structure  that stores the propterties of the Problem



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
    cost = 0;   % running Cost
    Cost_total = cost;
    
    % input of the previous time step
    u_last = [0;0];
    
    % the NMPC iteration FROM THE Book
    mpciter = 0;
    while(mpciter < mpciterations)
        % Step 1:
        % Obtain new initial value
        x0 = xmeasure;
        % Obtain matrices
        [Ad, Bd, x0fcT] = computeCurvedBicycleModel(x0, T, 0, Properties_Obj.l_f, Properties_Obj.l_r);
        myMatrices.A=Ad;
        myMatrices.B= Bd;
        myMatrices.x0fcT = x0fcT;
        
        % check if an accident can happen during the next predicted Horizon
        accident = zeros(N+1,1);
        for c =1:(N+1)
            if (x0(1)+Properties_Obj.l_f+ Properties_Obj.safety_distance>= x_refTV(1,mpciter+1)-Properties_Obj.l_r)&& ( x0(1)-Properties_Obj.l_r<=x_refTV(1,mpciter+1)+Properties_Obj.l_f+Properties_Obj.safety_distance) 
                accident(c) = 1;
            end
        end
 
        % Step 2 :
        %   Solve the optimal control problem
        t_Start = tic;
        [u_new, V_current, exitflag, output] = solveOptimalControlProblem ...
            (runningcosts, constraints, ...
             linearconstraints, system, ...
            N, T, x0, u0, options, x_ref(:,mpciter+1:mpciter+N+1),x_refTV, ...
        Properties_Obj , myMatrices, u_last,Determine_X_new, accident);

        t_Elapsed = toc( t_Start );
        
        % compute the predicted trajectory of the EV
        x_pred = dynamics(system,N,T,x0,u_new, Properties_Obj, myMatrices,Determine_X_new);
        
        %   Print solution
        printSolution(printHeader, printClosedloopData, ...
                      plotTrajectories, mpciter, x0, u_new, ...
                      exitflag, output, t_Elapsed, Properties_Obj, x_ref,x_pred, x_refTV(:,mpciter+1));


        %   Store closed loop data
        x = [ x xmeasure ]; 
        u = [ u u_new(:,1)];
        
        % Step 3:
        %  Apply control to process
        xmeasure = applyControl(system, T, x0, u_new, myMatrices);
        
        % add previous Cost value
        cost = cost + runningcosts(xmeasure, [u_last,u_new(:,1)], x_ref(:,1), Properties_Obj);
        Cost_total = [Cost_total cost];
     
        % Cost Plot
%         plotCost(Cost_total, mpciter);
% 
%         % Velocity Graph
%         plotVelocity(x, Properties_Obj, mpciterations, mpciter);
% 
%         % Position d
%         plot_d(mpciter, mpciterations,N, x);

        % plot Cars allure
        plotTrajectories( x, x0, x_ref,Properties_Obj,x_refTV(:,mpciter+1))
        
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
    options, x_ref,x_refTV, Properties_Obj, myMatrices, u_last,Determine_X_new, accident)
    % Set control and linear bounds
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    for k=1:N
        [Anew, bnew, Aeqnew, beqnew, lbnew, ubnew] = ...
               linearconstraints();
        A = blkdiag(A,Anew);
        b = [b, bnew];
        Aeq = blkdiag(Aeq,Aeqnew);
        beq = [beq, beqnew];
        lb = [lb, lbnew];
        ub = [ub, ubnew];
    end

    % Solve optimization problem
    [u, V, exitflag, output] = fmincon(@(u) costfunction(runningcosts, ...
    system, N,T, x0, ...
    u,x_ref,Properties_Obj, myMatrices, u_last,Determine_X_new), u0, A, b, Aeq, beq, lb, ...
    ub, @(u) nonlinearconstraints(constraints, ...
    system, N,T, x0, u, Properties_Obj,x_ref,x_refTV,myMatrices,Determine_X_new, accident), options);
        

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
                    N,T, x0, u, x_ref, Properties_Obj, ...
                    myMatrices, u_last,Determine_X_new)
    cost = 0;
    n = length(x0);
    x = zeros(n, N+1);
    % x: openloopsolution
    x = dynamics(system,N,T,x0,u, Properties_Obj, myMatrices,Determine_X_new); 
    cost = cost+ runningcosts(x(:,2), [u_last,u(:,1)], x_ref(:,2), Properties_Obj);
    
    for k=2:N
        cost = cost+runningcosts(x(:,k+1), u(:,k-1:k), x_ref(:,k+1), Properties_Obj);
    end
end

function [c,ceq] = nonlinearconstraints(constraints, ...
    system, ...
    N,T, x0, u, Properties_Obj,x_ref,x_refTV,myMatrices, Determine_X_new,accident)
    x = zeros(N+1, length(x0));

    %OpenloopSolution
    x = dynamics(system,N,T,x0,u, Properties_Obj, myMatrices, Determine_X_new);
    c = [];
    ceq = [];
    for k=2:N
        [cnew, ceqnew] = constraints(x(:,k),x_ref(:,k),x_refTV, k, Properties_Obj,accident);
        c = [c cnew];
        ceq = [ceq ceqnew];
    end
    % terminal constraints
    [cnew, ceqnew] = constraints(x(:,N+1),x_ref(:,N+1),x_refTV, N+1, Properties_Obj,accident);
    c = [c cnew];
    ceq = [ceq ceqnew];
end


%%%%%% OUTPUT %%%%%%%%%%%%%
function printSolution(printHeader, printClosedloopData, ...
             plotTrajectories, mpciter, x0, u, exitflag, output, t_Elapsed, ...
             Properties_Obj, x_ref, x , x_refTV)
    % Print results
    if (mpciter == 0)
        printHeader();
    end
    printClosedloopData(mpciter, u, x0, t_Elapsed);
    
    plotTrajectories(x, x0, x_ref,Properties_Obj,x_refTV ) 
end


function X = dynamics(system, N,T, x0, u0, Properties_Obj, myMatrices,Determine_X_new)
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
function plotVelocity(x, Properties_Obj, mpciterations, mpciter)
    figure(3);
    plot(0:mpciter,x(4,:),'r');
    axis([0 mpciterations 0 Properties_Obj.v_max]);
    title('Velocity  NMPC TEST ');
    xlabel('time (k)');
    ylabel('velocity v (m/s)');
end

function plotCost(Cost_total, mpciter)
   figure(4)
   plot(0:mpciter,Cost_total(1:mpciter+1),'r');
   title('Standard MPC');
   xlabel('time instant k');
   ylabel('Cost function value');
end

function plot_d(mpciter, mpciterations,N, x)
    figure(5);
    plot(0:mpciter,x(2,:),'r');
    axis([0 mpciterations+N+1 -5 5]);
    title('d  NMPC TEST ');
    xlabel('time (k)');
    ylabel('d (m)');
end
