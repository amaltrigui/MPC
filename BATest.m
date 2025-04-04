function Main
%% General
    clc;
    close all;
    clear all;
    addpath myfunctionsfolder;
    
    figure(1);
        title('Test Two Cars on the same Lane ');
        xlabel(' x in meter');
        ylabel(' y in meter');
        grid on;
        hold on;
        
        
    mpciterations = 50;     % the Length of the Simulation 
    N             = 12 ;    % prediction horizon  
    T             = 0.2;    % Sampling interval
    

%% EV and Lane in EV_system
    
    % Covariance Matrix, Gaussian Noise, Error Matrix
    covar_matrix = [ 0.15 0 ; 0 0.03];
    TV_system.covar_matrix = covar_matrix; 
    w = [0;0];
    w(1) = normrnd(0,covar_matrix(1,1));
    w(2) = normrnd(0,covar_matrix(2,2));
    TV_system.noise = w;
    
    % Beta = Risk Parameter
    TV_system.beta = 0.8;
    TV_system.extraSafetyDistance = 4;

    % Car Properties
    EV_system.l_r = 2;    
    EV_system.l_f = 2;  
    
    EV_system.Length = 5;
    EV_system.Width = 2;
    
    EV_system.v_max = 20; 

    % Matrices for(Cost function) from the Paper (see Page 6/7)
    EV_system.Q = diag([0 1 1 1]);
    EV_system.R = diag([0.33 5]);
    EV_system.S = diag([0.33 15]);

    EV_system.v_max = 20; 

    % Road Properties
    EV_system.Lane = [-3,3]; % Since LaneWidth is 3, the first lane is
    % defined between -3,0 and the second Lane is between 0,3

    EV_system.maxDeviation = 4; 
    EV_system.safety_distance = 2;

%% Initial values
    % Initial state EV X: (4*1)Vector (s,d,phi,v)
    
    % EV:
    % Velocity of the EV is 7m/S
    xmeasure      = [0.0; -1.5; 0.0; 9];  
    
    % Input:
    u0            = zeros(2,N);  
    
    % Initial state TV X: (4*1)Vector (x,vx,y,vy)

    % TV:
    % Velocity of Car 2 is 5m/S, smaller than velocity of EV --> to get the
    % Collision Case. TV is already at the Beginning of the Simulation 16 m ahead of EV.
    x_TV_measure  = [16; 5; -1.5; 0];  
       
%% reference trajectory

    % EV:
    % How is it wished, so that the car goes if no constaints are added 
    %              --> d 'lateral position' would stay 0(No deviation)
    %              --> phi 'Orientation' is also 0
    %              --> v would stay kostant 7 m/s
    %              --> s 'longitudinal position' goes from 0:v*t and t is 1 second for exemple ? Does
    %              it matter since Q(0) =0000--> No
    
    d           =   xmeasure(2)*ones(1,mpciterations+N+1);
    phi         =   zeros(1,mpciterations+N+1);
    v           =   xmeasure(4)*ones(1,mpciterations+N+1);
    x_ref = [0:mpciterations+N;d;phi;v];
    % TV:
    % Reference Trajectory of TV
    % assuming that the second car maintains the same velocity and same
    % lane as Reference Trajectory
    
        
     yTV         =   x_TV_measure(3)*ones(1,mpciterations+N+1);
     vyTV        =   zeros(1,mpciterations+N+1);
     vxTV        =   x_TV_measure(2)*ones(1,mpciterations+N+1);
     xTV         =   (0:1:mpciterations+N)*T*x_TV_measure(2)+ x_TV_measure(1);
     x_TV_ref    =   [xTV;vxTV;yTV;vyTV];   
    
     
%% offline Commands

    TV_system.sigma_e = cov_propagation(N, T, TV_system);
    var_s = [];
    e_s = [];
    var_d = [];
    e_d = [];
    
    for k = 1:N
        % s : longitudinal position cooresponds to x coordinate
        var_s_k = sqrt(TV_system.sigma_e(1,1,k));
        var_s = [var_s var_s_k] ;
        e_sk_k =  var_s_k*sqrt(-2*log(1-TV_system.beta));
        e_s = [e_s e_sk_k];
        
        % d : lateral position cooresponds to y coordinate
        var_d_k = sqrt(TV_system.sigma_e(3,3,k));
        var_d = [var_d var_d_k] ;
        e_dk_k =  var_d_k*sqrt(-2*log(1-TV_system.beta));
        e_d = [e_d e_dk_k];
    end
    TV_system.var_s = var_s;
    TV_system.e_s = e_s;
    
    TV_system.var_d = var_d;
    TV_system.e_d = e_d;
    TV_system.b_r = EV_system.Width + TV_system.extraSafetyDistance + TV_system.e_d ;
    TV_system.a_r_velocity_independent = EV_system.Length + TV_system.extraSafetyDistance + TV_system.e_s ;
    % it is missing the term a(X0, Xtv) that should be added after findind X0
    TV_system.a_r_velocity_dependent = TV_system.a_r_velocity_independent ;

%% Optimization
    nmpc(@runningcosts, @constraints, ...
          @linearconstraints, @system, ...
         mpciterations, N, T, xmeasure, u0, x_ref, x_TV_measure, x_TV_ref, @Determine_X_new, ...
            EV_system,TV_system, @computeCurvedBicycleModel, @printHeader, @printClosedloopData, @plotTrajectories, @TV_prediction, @TV_dynamics);

    rmpath('./myfunctionsfolder');

end


function cost = runningcosts(x, u, x_ref,EV_system)
%% From the Paper , see Formel (25)
   u_prev = u(:,1);
   u_curr = u(:,2);
   
   %% .' is the right Syntax
   cost = (x-x_ref).'*EV_system.Q*(x-x_ref) + u_curr.'*EV_system.R*(u_curr) + (u_curr-u_prev).'*EV_system.S*(u_curr-u_prev);
end


function [c,ceq] = constraints(x, x_ref, x_TV, k, EV_system, TV_system)
%% constraints
    

    % the Rectangular TV safety area is defined by :
    b_r_k = TV_system.b_r(k);
    a_r_k = TV_system.a_r_velocity_dependent(k);


    % velocity of EV positive and smaller than Vmax
    vmax = x(4)- EV_system.v_max ;  % v<vmax
    v_positive = -x(4);                  % v>0
    
    % do not deviate so much from reference Trajectory 
    maxDev = EV_system.maxDeviation; 
    R = abs(x(2)-x_ref(2))-maxDev;
    
    % stay in the Road/ one of the two Lanes
    EV_Width =  EV_system.Width ;
    upper_limit = EV_system.Lane(2); 
    lower_limit = EV_system.Lane(1); 
    
    l1 = -x(2) + EV_Width/2 + lower_limit;
    l2 = x(2) + EV_Width/2 - upper_limit; 
    % last Element of c is -2, which is always negative but in case of
    % accident detected it changes to the constraint that prevent it
    c   = [vmax; v_positive ; R ; l1; l2; -2];

    
    % Safety Constraint Formel 11
    safetyDistance = (EV_system.Length/2) + TV_system.e_s(1,k) + TV_system.extraSafetyDistance;
    
    new_SMPC_Constarint = constraint_Genereation (x, x_TV, EV_system, safetyDistance, a_r_k, b_r_k);
    % new way: 
    c(end) = new_SMPC_Constarint;
    
    %% old way
    % check if an accident can happen at step k (out of N horizon steps)
%     accident_boolean = (x(1)+EV_system.l_f+ EV_system.safety_distance>= x_TV(1)-EV_system.l_r) && ( x(1)-EV_system.l_r<=x_TV(1)+EV_system.l_f+EV_system.safety_distance);
%     if(accident_boolean)
%         A = x_TV(3)+ EV_Width- x(2);
%         c(end)=A;
%     end
          
    ceq = [];
end

function new_SMPC_Constraint = constraint_Genereation(x, x_TV, EV_system, safetyDistance, a_r_k, b_r_k)
    delta_x = x(1)-x_TV(1);
    r_close = 90;
    
    if abs(delta_x) >= safetyDistance
        new_SMPC_Constraint = -1; % No additional constraint, EV is safe and far away from TV
    elseif (delta_x < 0 && - delta_x >r_close )
        c_2_x_k = x_TV(1)- a_r_k/2;
        new_SMPC_Constraint = x(1) - c_2_x_k; % stay behind the TV
    elseif (delta_x > 0 && delta_x > r_close)
        c_1_x_k = x_TV(1) + a_r_k/2;
        new_SMPC_Constraint = - x(1) + c_1_x_k; % stay in front of the TV
    elseif ( (x(4) == x_TV(3)) && (- delta_x <= r_close) )
        % accident fall : EV is so close and behind to TV and in the same lane
        c_4_y_ev = x(2)- EV_system.Width/2;
        c_2_y_k =  x_TV(3) + b_r_k/2;
        c_4_x_ev = x(1) +  EV_system.Length/2;
        c_2_x_k =  x_TV(1) - a_r_k/2;
        % Constraint that will push EV to the left lane
        new_SMPC_Constraint = -x(2) + max([0, (c_4_y_ev-c_2_y_k)/(c_4_x_ev-c_2_x_k)])* x(1) + c_2_y_k- x(1)*c_4_x_ev;
    elseif ((x(4) >= x_TV(3)) && (abs(delta_x) < r_close))
        % EV is on the right to TV but x-coordinates are close --> goes
        % parallel to TV and thus should maintain his Lane and not deviate
        % back to right lane
        c_2_y_k = x_TV(3) + b_r_k/2;
        new_SMPC_Constraint = -x(2) + c_2_y_k;
        
    end

end

function [A, b, Aeq, beq, lb, ub] = linearconstraints()
%% Input Constraints
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = [-9; -0.52];   % Min input bound (paper)
    ub  = [5;0.52];      % Max input bound
end



function x_new = Determine_X_new(x, u, EV_system)
 %% formel (18) Paper
 
 s = x(1); d = x(2); phi = x(3); v = x(4); a = u(1); delta = u(2);
  
 l_r = EV_system.l_r;  l_f = EV_system.l_f;
 
 alpha = atan((l_r*tan(delta))/(l_f + l_r));
 s_new = v*cos(phi + alpha);
 d_new = v*sin(phi + alpha);
 phi_new = (v*sin(alpha))/l_r;
 v_new = a;

 x_new = [s_new;d_new;phi_new;v_new];

end


 function y = system(x, u, T, x_curr, myMatrices)
 %% system Equation
    A_d = myMatrices.A;
    B_d = myMatrices.B;
    x0fcT = myMatrices.x0fcT;

    y = x0fcT + A_d*(x - x_curr) + B_d*u;
 end
 

 %% print results
 function printHeader()
    fprintf('   k  |      a(k)        delta(k)        s        d        phi        v     Time\n');
    fprintf('--------------------------------------------------------------------------------\n');
 end

         
function printClosedloopData(mpciter, u, x, t_Elapsed, x_TV)
%% from Web
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+6.3f \n', ...
             mpciter, u(1,1),u(2,1), x(1), x(2),x(3),x(4), t_Elapsed);
    
    fprintf('\n ');     
%     fprintf(' Values of TV state : \n ') ;
%     x_TV
end


function plotTrajectories( x, x0, x_ref, EV_system,x_TV)
 %% Visualize Cars   
    % Assuming that the Width and Length of EV is equal to TV

    L = EV_system.Length ;  %  length of car
    B = EV_system.Width  ;  %  width of car

    C = [0. 0. ; L 0. ; L B ; 0. B; 0 0.] ; % Car Center coordinates ( from Internet copied--> stil need to check it butit  works )

    borders_lanes = EV_system.Lane;
    middle_Road = (borders_lanes(1)+borders_lanes(2))/2; % the two lanes are limited


    % plot street
    plot(linspace(-2,90,100),borders_lanes(1)*ones(100),'b')  % upper Borderie
    plot(linspace(-2,90,100),borders_lanes(2)*ones(100),'b')  % lower Borderie
    plot(linspace(-2,90,100),middle_Road*ones(100), '--')    % middle of Road

    % plot x_Ref EV
    plot(linspace(x_ref(1,1),x_ref(1,length(x_ref)),100),linspace(x_ref(2,1),x_ref(2,length(x_ref)), 100),'r')
    
    % plot EV trajectory
    h8= plot(C(:,1)+x0(1)-L/2,C(:,2)+x0(2)-B/2,'black');

    % plot TV trajectory
    h9 =plot(C(:,1)+x_TV(1)-L/2,C(:,2)+x_TV(3)-B/2,'green');

    axis([-2 90 -25 25]); 
    pause(0.3)
    set(h8,'Visible','off')
    set(h9,'Visible','off')

end


function [Ad, Bd, x0fcT] = computeCurvedBicycleModel(x, T, k, lf, lr)
%% CurvedBicycleModel Paper/Tommaso

    d = x(2);
    v = x(4);
    c0 = cos(x(3));
    s0 = sin(x(3));
    alpha1 = lr/(lr+lf);
    
    w0 = 1-k*d;
    w1 = 1/w0;
    lambda = v*k*w1;

    z0 = sqrt(1-5*c0^2);  %   z1-z2=2*z0
    z1 = s0+z0;
    z2 = s0-z0;
    e1 = exp(lambda*z1*T/2);
    e2 = exp(lambda*z2*T/2);
    z3 = e1-e2;
    z4 = z1*e2-z2*e1;
    z5 = z1*e1-z2*e2;

    if abs(z0)<1e-5
        z6 = lambda*T;
        z7 = 2;
        z8 = 2*(1+lambda*s0*T);
        z9 = lambda*T^2/2;
        z10 = 2*T;
        z11 = 2*T+lambda*s0*T^2;
    else
        z6 = z3/z0;
        z7 = z4/z0;
        z8 = z5/z0;
        if lambda==0
            z9 = 0;
            z10 = 2*T;
            z11 = 2*T+lambda*s0*T^2;
        else
            if z1*z2==0
                if z1==0
                    z9 = (T-2*(e2-1)/lambda/z2)/z0;
                    z10 = -z2/z0*T;
                else
                    z9 = (2*(e1-1)/lambda/z1-T)/z0;
                    z10 = z1/z0*T;
                end
            else
                z9 = 2*((e1-1)/z1-(e2-1)/z2)/z0/lambda;
                z10 = 2*((e2-1)*z1/z2-(e1-1)*z2/z1)/z0/lambda;
            end
            z11 = 2*z6/lambda;
        end
    end

    if lambda==0
        a13 = -v*s0*T;
        a14 = T*c0*w1;
        a23 = v*c0*T;
        a24 = T*s0;
        a34 = -k*c0*T*w1;
        b11 = c0*T^2/2*w1;
        b21 = s0*T^2/2;
        b31 = -k*T^2*c0/2*w1;
        b12 = -(1+v*T/2/lr)*v*alpha1*s0*T*w1;
        b22 = (1+v*T/2/lr)*v*T*c0*alpha1;
        b32 = v*alpha1*T/lr;
    else
        if c0==0
            a14 = 0;
            a24 = s0*T;
            a34 = 0;
            b11 = 0;
            b21 = s0*T^2/2;
            b31 = 0;
        else
            a14 = (s0*(1-z8/2)+z6)/(v*k*c0);
            a24 = (-1+s0*c0^2*z6+z7/2)/lambda/c0^2;
            a34 = (s0*(z8/2-1)-z6)/v/c0;
            b11 = (s0*(T-z11/2)+z9)/(v*k*c0);
            b21 = (-T+s0*c0^2*z9+z10/2)/lambda/c0^2;
            b31 = (s0*(z11/2-T)-z9)/v/c0;
        end
        a13 = (1-z8/2)/k;
        a23 = w0*c0*z6/k;
        w2 = (w0/k/lr+s0);
        b12 = (-T*s0+c0^2*z9+(T-z11/2)*w2)*v*alpha1*w1;
        b22 = (z10/2+z9*w2)*v*c0*alpha1;
        b32 = (z11*w2/2-z9*c0^2)*lambda*alpha1;
    end

    Ad = [1,    c0*w1*z6,    a13,         a14;
          0,    z7/2,             a23,         a24;
          0,    -k*c0*z6*w1, z8/2,   a34;
          0,    0,                    0,          1];

    Bd = [b11,   b12;
          b21,   b22;
          b31,   b32;
          T,    0];
      
    q1 = 1/(1-k*x(2));    %   1/(1-k*d)
    q2 = cos(x(3));        %   cos(phi)
    q3 = x(4)*q2*q1;       %   s_dot = v*cos(phi)/(1-k*d)
    fcq = [q3; x(4)*sin(x(3)); -k*q3; 0];
    x0fcT = x+T*fcq;
end


%%%%%%%% TV CAR %%%%%%%

% for prediction no noise
% for real trajectory with noise


% xTV = [x,vx,y,vy] 4 x 1 
%--> see paper ,copy values

function [x_TV_withnoise, x_TV] = TV_dynamics(xTVk, x_TV_ref_k, T, w) 
%% x_TV_ref_k is the k-st eintrag in x_TV_ref
    % system matrices
    A = [1 T 0 0 ; 0 1 0 0; 0 0 1 T; 0 0 0 1];
    t = (T^2)/2;
    B = [t 0; T 0; 0 t; 0 T];
    
    % stabilizing feedback matrix
    K = [ 0 -0.55 0 0 ; 0 0 -0.63 -1.15];
    
    % Prediction of TV without noise
    x_TV = (A+B*K)*xTVk - B*K*x_TV_ref_k;
    
    % Add Gaussian Noise
    x_TV_withnoise = x_TV + B*w ;
end

function X_TV = TV_prediction(x_TV_0, x_TV_ref, N, T, w)
%%  predict X_TV  for N Horizon
    X_TV(:,1) = x_TV_0;
    for k=1:N
       [~,xtv]  = TV_dynamics(X_TV(:,k),x_TV_ref(:,k), T,w);  
        X_TV(:,k+1)= xtv;   
    end
end

function sigma_e = cov_propagation(N, T, TV_system)
%%
    % covariance matrix propagation

    w_cov = TV_system.covar_matrix;
   
    %define them first , store them in Object SystemTV, pass it als
    %parameter
    
    A = [1 T 0 0 ; 0 1 0 0; 0 0 1 T; 0 0 0 1];
    t = (T^2)/2;
    B = [t 0; T 0; 0 t; 0 T]; 
    K = [ 0 -0.55 0 0 ; 0 0 -0.63 -1.15];

    r = A+B*K;

    % initial covariance matrix
    sigma_e = zeros(4,4,N);

    % update covariance matrix for each step
    for i = 2:N
        sigma_e(:,:,i) = r*sigma_e(:,:,i-1)*r' + B*w_cov*B';
    end
    
end

