function [Ad, Bd, x0fcT] = computeCurvedBicycleModel(x, T, k, lf, lr)

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
