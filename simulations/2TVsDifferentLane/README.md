**Description**

Second Test : 
* Road is composed of two lanes (only one direction) of width 3 m (each one)
* EV starts from Start (0 m, -1.5 m == Center of the right Lane)
* 2 TVs: the first TV starts from  (40 m, -1.5 m ), 50 m in front of EV but has a smaller velocity (7 m/s) than the EV Velocity (15 m/s)
          the second TV starts from (50 , 1.5 m), is in the left Lane  from EV but has a smaller velocity (10 m/s) than the EV Velocity (15 m/s).
          The second TV represents an Obstacle in front of EV by overtaking the first TV, Thus EV will slow down on the right Lane until the safety Distance between them is respected, then overtakes it to the left
* At a certain time, since EV is faster than TV, EV will overtake the first TV from the left than returns to the right Lane  (The right Lane is always the reference Lane in the Road)
* Plots :

Velocity           |  Trajectory 
:-------------------------:|:-------------------------:
![](https://github.com/amaltrigui/MPC/blob/241a3506b683e53d4b8ca3eaa003fa2b22664051/simulations/2TVsDifferentLane/velocity.jpg )  |  ![](https://github.com/amaltrigui/MPC/blob/241a3506b683e53d4b8ca3eaa003fa2b22664051/simulations/2TVsDifferentLane/d.jpg)

current Cost           |  Total Cost 
:-------------------------:|:-------------------------:
![](https://github.com/amaltrigui/MPC/blob/241a3506b683e53d4b8ca3eaa003fa2b22664051/simulations/2TVsDifferentLane/cost.jpg )  |  ![](https://github.com/amaltrigui/MPC/blob/241a3506b683e53d4b8ca3eaa003fa2b22664051/simulations/2TVsDifferentLane/totalCost.jpg)



