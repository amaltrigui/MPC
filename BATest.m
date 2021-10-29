function Main
%% General
    clc;
    close all;
    clear all;
    addpath myfunctionsfolder;
    
%     figure(1);
%         title('Test One Car ');
%         xlabel(' x in meter');
%         ylabel(' y in meter');
%         grid on;
%         hold on;
        
        
    mpciterations = 50;     % the Length of the Simulation 
    N             = 6;      % prediction horizon  
    T             = 0.2;    % Sampling interval
    

%% EV and Lane Properties_Objeters
    Properties_Obj.l_r = 2;    
    Properties_Obj.l_f = 2;    
    
    % Matrices for(Cost function) from the Paper (see Page 6/7)
    Properties_Obj.Q = diag([0 1 1 1]);
    Properties_Obj.R = diag([0.33 5]);
    Properties_Obj.S = diag([0.33 15]);

    Properties_Obj.v_max = 13; 

    Properties_Obj.Lane = [-2,2]; % lane width is 4
 

%% Initial values
    % Initial state X: (4*1)Vector (s,d,phi,v)
    % Car in Q(0,0) in Coorinate-System
    % Car has 0 rad as orientation of Veh bzgl Road
    % Car has an initial Velocity of 3 m/s
    xmeasure      = [0.0; 0.0; 0.0; 3];  
    
    % Initial Input
    % u0 : (2*1)
    % accelearion is 0 m/s^2
    % steering angle is 0 rad  --> Car has a konstant velocity
    u0            = zeros(2,N);  
    
%% reference trajectory
    % How would the car goes if no constaints are added --> No acceleration
    % or deviation --> d 'lateral position' would stay 0(No deviation)
    %              --> phi 'Orientation' is also 0
    %              --> v would stay kostant 9 m/s
    %              --> s 'longitudinal position' goes from 0:v*t and t is 1 second for exemple ? Does
    %              it matter since Q(0) =0000
    d = zeros(1,mpciterations+N+1);
    phi =  zeros(1,mpciterations+N+1);
    v = 3*ones(1,mpciterations+N+1);
    v(10:29)=7;  % change to reference Velocity between Time 10 and 27 (instant) to 7
    x_ref = [0:mpciterations+N;d;phi;v];
    
%% Optimization

    Properties_Obj.steps = 1;
    nmpc(@runningcosts, @constraints, ...
          @linearconstraints, @system, ...
         mpciterations, N, T, xmeasure, u0, x_ref,@Determine_X_new, ...
            Properties_Obj, @computeCurvedBicycleModel, @printHeader, @printClosedloopData, @plotTrajectories);

    rmpath('./functions');

end


function cost = runningcosts(x, u, x_ref,Properties_Obj)
%% Paper Formel (25)

   u_prev = u(:,1);
   u_curr = u(:,2);
  
   cost = (x-x_ref).'*Properties_Obj.Q*(x-x_ref) + u_curr.'*Properties_Obj.R*(u_curr) + (u_curr-u_prev).'*Properties_Obj.S*(u_curr-u_prev);
end


function [c,ceq] = constraints(x, x_ref, k, Properties_Obj)
%% constraints
    
    vmax = x(4)- Properties_Obj.v_max;  % v<vmax
    v_positive = -x(4);        % v>0
    c   = [vmax;v_positive];
    ceq = [];
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



function x_new = Determine_X_new(x, u, Properties_Obj)
 %% Dynamics of the non linearized system: update current state but not for prediction
 % formel (18) Paper
 
 s = x(1); d = x(2); phi = x(3); v = x(4); a = u(1); delta = u(2);
  
 l_r = Properties_Obj.l_r;  l_f = Properties_Obj.l_f;
 
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

function printClosedloopData(mpciter, u, x, t_Elapsed)
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+6.3f', ...
             mpciter, u(1,1),u(2,1), x(1), x(2),x(3),x(4), t_Elapsed);
end

function plotTrajectories( x, x0, Properties_Obj, x_ref)      
% plot road
    l_r = Properties_Obj.l_r;
    l_f = Properties_Obj.l_f;
    L = l_r + l_f;  %  length of car
    B = L/2;        %  width of car
    
    C = [0. 0. ; L 0. ; L B ; 0. B; 0 0.] ; % Car Center coordinates (Internet copied)
    
    
    y_lim = Properties_Obj.Lane;
    y_middle = (y_lim(1)+y_lim(2))/2; % lane cut in 2
        
%     plot(linspace(-2,40,100),y_lim(1)*ones(100),'b')  % upper Borderie
%     plot(linspace(-2,40,100),y_lim(2)*ones(100),'b')  % lower Borderie
%     plot(linspace(-2,40,100),y_middle*ones(100), '--')    % middle of Lane

% plot x_Ref
%     plot(linspace(x_ref(1,1),x_ref(1,length(x_ref)),100),linspace(x_ref(2,1),x_ref(2,length(x_ref)), 100),'r')
% plot car
%     h4= plot(C(:,1)+x0(1)-L/2,C(:,2)+x0(2)-B/2,'black')

% adjust axis
%     axis([-2 40 -5 7]);
%     pause(1)
%     set(h4,'Visible','off')
    
end


%% CurvedBicycleModel Paper/Tomasso
function [Ad, Bd, x0fcT] = computeCurvedBicycleModel(x, T, k, lf, lr)

    x
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







