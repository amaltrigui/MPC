**Description**

Second Test : 
* Road is composed of two lanes (only one direction) of width 3 m (each one)
* EV starts from Start (0 m, -1.5 m == Center of the right Lane)
* 2 TVs: the first TV starts from  (50 m, -1.5 m ), 50 m in front of EV but has a smaller velocity (7 m/s) than the EV Velocity (15 m/s)
          the second TV starts from (65m , -1.5 m), 65 m in front of EV but has a smaller velocity (7 m/s) than the EV Velocity (15 m/s)
* At a certain time, since EV is faster than TV, EV will overtake the 2 TVs from the left than returns to the right Lane (The right Lane is always the reference Lane in the Road)
* Plots :

Velocity           |  Trajectory 
:-------------------------:|:-------------------------:
![](https://github.com/amaltrigui/MPC/blob/f16b651b305032b727072d7b9fbb1b7869ba38ca/simulations/SecondTest2TVs/velocity.jpg )  |  ![](https://github.com/amaltrigui/MPC/blob/f16b651b305032b727072d7b9fbb1b7869ba38ca/simulations/SecondTest2TVs/d.jpg)

current Cost           |  Total Cost 
:-------------------------:|:-------------------------:
![](https://github.com/amaltrigui/MPC/blob/f16b651b305032b727072d7b9fbb1b7869ba38ca/simulations/SecondTest2TVs/cost.jpg )  |  ![](https://github.com/amaltrigui/MPC/blob/f16b651b305032b727072d7b9fbb1b7869ba38ca/simulations/SecondTest2TVs/totalCost.jpg)



