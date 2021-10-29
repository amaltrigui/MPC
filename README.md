# MPC

shared Code for MPC

#### step 1: non linear MPC, simple Simulation:
As **Parameters** I chose:
* mpciterations = 50;     
* N             = 6;      
* T             = 0.2; 
* X_reference = [s_ref, d_ref, phi_ref, v_ref], where:
    - s_ref : random for the moment 
    - d_ref : 0
    - phi : 0
    - v_ref : for k 10:29 (out of mpciterations+N+1 )is v = 7m/s otherwiese 3 m/s =v_0 --> Just to check the Velocity Curve

* X_0 = [0,0,0,3]

The plot of Velocity looks like:

!<img src="https://github.com/amaltrigui/MPC/blob/main/Plots/v1.PNG" width="400" height="400" />
