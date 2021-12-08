**Description**

First Test : 
* Road is composed of two lanes (only one direction) of width 3 m (each one)
* EV starts from Start (0 m, -1.5 m == Center of the right Lane)
* TV starts from Start (60 m, -1.5 m ), 60 m in front of EV but has a smaller velocity (10 m/s) than the EV Velocity (20 m/s)
* At a certain time, since EV is faster than TV, EV will overtake TV from the left than returns to the right Lane (The right Lane is always the reference Lane in the Road)
* Plots :

Velocity           |  Trajectory 
:-------------------------:|:-------------------------:
![](https://github.com/amaltrigui/MPC/blob/dbaac28238c5cca27455e6c86f9402804b258fb2/simulations/SimpleTest1TV/velocity.jpg )  |  ![](https://github.com/amaltrigui/MPC/blob/ef3cc8539deb65241de16c830126b030d448b9b0/simulations/SimpleTest1TV/d.jpg)

current Cost           |  Total Cost 
:-------------------------:|:-------------------------:
![](https://github.com/amaltrigui/MPC/blob/ef3cc8539deb65241de16c830126b030d448b9b0/simulations/SimpleTest1TV/cost.jpg )  |  ![](https://github.com/amaltrigui/MPC/blob/ef3cc8539deb65241de16c830126b030d448b9b0/simulations/SimpleTest1TV/totalCost.jpg)


