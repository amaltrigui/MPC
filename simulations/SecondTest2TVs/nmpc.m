function [x, u, totalcost] = nmpc(runningcosts, ...
              constraints, ...
              linearconstraints, system, ...
              mpciterations, N, T,  xmeasure, u0, x_ref, xTVs, ...
              Determine_X_new, EV_system,TV_system, computeCurvedBicycleModel, ...
              printHeader, printClosedloopData, plotTrajectories ,TV_prediction, TV_dynamics, constraint_Coeff_Genereation, writerObj, RefrenceTV)

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
    x = [];  % real Position of the EV
    u = [];  % Inputs
    u_last = [0;0]; % input of the previous time step

    totalcost = 0;
    cost_array = 0;
    cost = 0;   % the running Cost
    Cost_total = cost; % total Cost of the whole Simulation
    
    x_TV_predicted = []; % real position of TV
    w = TV_system.noise;
        
    % the NMPC iteration (FROM THE Book)
    mpciter = 0;
    while(mpciter < mpciterations)
        
        %%%%%%%%%%%%%%%%%%%%%%%% Step 1: Get new/current Values %%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Get new current value of x ( first x value of the N-Horizon X)
        x0 = xmeasure;       
        
        % Get matrices from Curved Bicycle Model Function
        [Ad, Bd, x0fcT] = computeCurvedBicycleModel(x0, T, 0, EV_system.l_f, EV_system.l_r);
        myMatrices.A=Ad;
        myMatrices.B= Bd;
        myMatrices.x0fcT = x0fcT;
        
        % Predict Trajectory of TVs  for N horizon
        
        % at the very beginning, (mpciter==0), x_TV_measure is the first
        % actual state of TV (right at the beginning of simulation). After
        % each iteration, x_TV_measure will be updated to the new first
        % element of the predicted N.size Vector X_TV    
        
        xTVcurr = xTVs;
        % for loop 
        for j = 1:length(xTVs(1,:))
            xTV0 = xTVs(:,j); % 4*1
            xTVref = RefrenceTV(xTV0, T, N); % 4* (N*1)
            X_TV(:,:,j) = TV_prediction(xTV0,xTVref,N,T,w); % 4* (N*1) * j
            
            [x_TV_measure_withnoise,x_TV_measure] = TV_dynamics(xTV0, xTVref(:,1),T, w);
            xTVs(:,j) = x_TV_measure;
        end
        
      
        % START OPC Time
        t_Start = tic;
        [u_new, V_current, exitflag, output] = solveOptimalControlProblem ...
            (runningcosts, constraints, ...
             linearconstraints, system, ...
             N, T, x0, u0, options, x_ref(:,mpciter+2:mpciter+N+1),X_TV(:,2:N+1,:),xTVcurr, ...
        EV_system , TV_system, myMatrices, u_last,Determine_X_new, constraint_Coeff_Genereation);
        % END OPC Time
        t_Elapsed = toc( t_Start );
        
        % the current new predicted trajectory of the EV N+1 sized -->
        % needed in plot
        x_pred = computeOpenloopSolution(system, N, T, x0, u_new, myMatrices);
        
        %  Print solution
        printSolution(printHeader, printClosedloopData, ...
                      plotTrajectories, mpciter, x0, u_new, ...
                      exitflag, output, t_Elapsed, EV_system, x_ref, x_pred, x_TV_measure);


        %  Store closed loop data
        x = [ x xmeasure ]; 
        u = [ u u_new(:,1)];
        
        %%%%%%%%%%%%%%%%%%%%%% Step 3: update trajectory and inputs %%%%%%%%%%%%%%%%%%%%
        %  Apply control --> obatin new x_measure
        xmeasure = applyControl(system, T, x0, u_new, myMatrices);
        
        % update previous Cost value
        run_cost = runningcosts(xmeasure, [u_last,u_new(:,1)], x_ref(:,1), EV_system);
        cost_array = [cost_array run_cost];
        
        cost = cost + run_cost;
        Cost_total = [Cost_total cost];
        
        totalcost =  cost;
        
        % plot graphs
          plot_graphs(Cost_total,cost_array, mpciter, mpciterations, EV_system,x,N);


        % plot Cars allure
%          plotTrajectories( x_pred, x0, x_ref,EV_system, TV_system, xTVs, X_TV, N, writerObj , mpciter == mpciterations-1)
        
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
    options, x_ref,x_TV,xTV0, EV_system,TV_system, myMatrices, u_last,Determine_X_new, constraint_Coeff_Genereation)
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
    
        % not so sure !!!
        a_r_k = (-1/a_min_longitudinal * max([0,(x0(4)^2- x_TV(2,k,1)^2)])); % /2
        a_r = [a_r a_r_k];
    end

    % update a_r 
    TV_system.a_r_velocity_dependent = TV_system.a_r_velocity_independent + a_r;

    % find coefficients of constraints qs, qd, qt
    ntv = length(xTV0(1,:));
    qx = zeros(ntv,N);
    qy = zeros(ntv,N);
    qt = zeros(ntv,N);
    
    for j = 1:ntv
        
        delta_S_X_0 = x0(1)- xTV0(1,j);  % delta xEV0, xTV0

        
        for k= 1:N
            [qx_k, qy_k, qt_k] = constraint_Coeff_Genereation(delta_S_X_0, x0(1),x0(2), x_TV(:,:,j), EV_system, TV_system, k);
             qx(j,k) =  qx_k;
             qy(j,k) =  qy_k;
             qt(j,k) =  qt_k;
        end
    end

    % Solve optimization problem
    [u, V, exitflag, output] = fmincon(@(u) costfunction(runningcosts, ...
    system, N,T, x0, ...
    u,x_ref,EV_system, myMatrices, u_last,Determine_X_new), u0, A, b, Aeq, beq, lb, ...
    ub, @(u) nonlinearconstraints(constraints, ...
    system, N,T, x0, u, EV_system, TV_system, x_ref,x_TV,myMatrices,Determine_X_new, qx, qy, qt), options);
        
%     keyboard();

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
        
    % x: openloopsolution
    x = computeOpenloopSolution(system, N, T, x0, u, myMatrices);
    cost = cost+ runningcosts(x(:,2), [u_last,u(:,1)], x_ref(:,1), EV_system);
    
    for k=2:N
        cost = cost+runningcosts(x(:,k+1), u(:,k-1:k), x_ref(:,k), EV_system);
    end
end

function [c,ceq] = nonlinearconstraints(constraints, ...
    system, ...
    N,T, x0, u, EV_system, TV_system, x_ref,x_TV,myMatrices, Determine_X_new, qx, qy, qt)

    x = zeros(length(x0),N+1);

    %OpenloopSolution
    x = computeOpenloopSolution(system, N, T, x0, u, myMatrices);
    c = [];
    ceq = [];
    
    for j =1:length(qx(:,1))
        for k=1:N
            [cnew, ceqnew] = constraints(x(:,k+1),x_ref(:,k), k, EV_system, TV_system, qx(j,k), qy(j,k), qt(j,k));
            c = [c cnew];
            ceq = [ceq ceqnew];
        end
    end
    % terminal constraints
%     [cnew, ceqnew] = constraints(x(:,N+1),x_ref(:,N+1),x_TV(:,N+1), N+1, EV_system, TV_system);
%     c = [c cnew];
%     ceq = [ceq ceqnew];
end


function printSolution(printHeader, printClosedloopData, ...
                       plotTrajectories, mpciter, x0, u,...
                       exitflag, output, t_Elapsed, EV_system, x_ref, x , x_TV)
%% Output   
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
end


%% Utils

function  plot_graphs(Cost_total,cost_array, mpciter, mpciterations, EV_system,x,N)

    Folder=pwd;
    [PathStr,FolderName]=fileparts(Folder);
    DataFolder=['plotsAndGraphs-',FolderName];
    if ~exist(DataFolder, 'dir')
        mkdir(DataFolder);
    end

    plotCost(Cost_total,cost_array, mpciter,DataFolder);
    plotVelocity(x, EV_system, mpciterations, mpciter,DataFolder);
    plot_d(mpciter, mpciterations,N, x,DataFolder);
    
end
function plotVelocity(x, EV_system, mpciterations, mpciter,DataFolder)
    figure(3);
    plot(0:mpciter,x(4,:),'r');
    axis([0 mpciterations 0 EV_system.v_max]);
    title('Velocity  NMPC TEST ');
    xlabel('time (k)');
    ylabel('velocity v (m/s)');
    saveas(figure(3),[DataFolder '/velocity.jpg']);

end

function plotCost(Cost_total,cost_array, mpciter,DataFolder)
   figure(4)
   plot(0:mpciter,Cost_total(1:mpciter+1),'r');
   title('Standard MPC');
   xlabel('time instant k');
   ylabel(' total Cost function value');
   saveas(figure(4),[DataFolder '/totalCost.jpg']);

   
    figure(9)
   plot(0:mpciter,cost_array(1:mpciter+1),'r');
   title('Standard MPC');
   xlabel('time instant k');
   ylabel('Cost function value');
   saveas(figure(9),[DataFolder '/cost.jpg']);

   
end

function plot_d(mpciter, mpciterations,N, x,DataFolder)
    figure(5);
    plot(0:mpciter,x(2,:),'r');
    axis([0 mpciterations+N+1 -5 5]);
    title('d  NMPC TEST ');
    xlabel('time (k)');
    ylabel('d (m)');
    
    
    saveas(figure(5),[DataFolder '/d.jpg']);
end


