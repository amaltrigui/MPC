# Simulations and results

shared Code for MPC

### First Example ) : 3 TVs + 1 EV:
#### High Level SMPC:
-> in this case, EV couldn t return back right because of TV3(green)

g

https://user-images.githubusercontent.com/50722211/151923123-5b159217-e41c-4058-8fd1-1d5e1e80657b.mp4


```diff
! Comparison cost =  382.5916
- total computation time: 62.8528 s
- total Low Level computation time : 56.5817 s
```

#### State Machine :
It depends on the chosen Parameters r_close (minimale distance between 2 vehicles that affects on the switching rules):

*1st : choose r_close = 60m ->  in this case the state machine behave almost similarly to Hierarchical MPC.
EV finds itself obliged to decrease the velocity after making the decision to deviate left and facing the second TV


https://user-images.githubusercontent.com/50722211/151923224-93421753-d639-4f2d-91e8-61bb948504a4.mp4


```diff
! Comparison cost =  575. 946
- total computation time: 39.5425 s
- Low Level compuutation time :39.2841 s
```
*2nd : choose r_close = 40m <1st r_close


```diff
! Comparison cost =  2567.6
- total computation time: 
- Low Level compuutation time :
```



### Second Example ) : 2 TVs + 1 EV:
#### High Level SMPC:
-Velocity_ev = 20m/s
-x_TV_measure  = [140; 7; -1.5; 0];  
-x_TV2_measure  = [150; 15; 1.5; 0];  15<20m/S
->Light decreasing in the velocity of EV because of TV2 when being in the left lane(from 20m/s to average 17m/s)

https://user-images.githubusercontent.com/50722211/151909788-3a5514be-f06e-44c8-8578-24078a49280c.mp4
#### State Machine :
video previously added (presentation)


### Simple Example :1TV + 1 EV
tried 5 different covariance matrices for the gaussian noise-->  the state machine maneuver develops a bigger comparison cost with more(increasing) noise than the one from HL SMPC.

Cost For the mpc: [4.5987    4.5987    4.5987    4.5987    4.5987]

Cost For the state machine: [ 15.9608   16.5978   17.1390   17.3863 18.1142]

![plotcosts](https://user-images.githubusercontent.com/50722211/151910280-99685dc1-271e-4900-be7a-e940ae0c5a1d.jpg)


### Run Simulation for n= 20 times
HL-MPC: check results (basically cost) each time for  a certain noise cov_matrice = [ 0.15 0 ; 0 0.03] and  mpciterations = 100 

--> static results

![untitled](https://user-images.githubusercontent.com/50722211/151920023-0a832c60-7400-47ef-920b-3aa7df4184e6.jpg)

State Machine: same check

total cost is highly increased 274.4603>>> 15.,,,

![untitledk](https://user-images.githubusercontent.com/50722211/151920840-c9bfd957-e00e-4faf-ada2-b510776c01df.jpg)

