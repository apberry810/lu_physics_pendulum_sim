%%%%%%-----------Forced-Damped-Pendulum-----%%%%%%%

These files were written and tested on Python 3.7.6
This is a simulation for a forced, damped, non-linear pendulum with arbitrary constant coefficients.
The pendulum has the equation of the form:
d^2(theta)/dt^2+B(d(theta)/dt+g*sin(theta)/l=A*cos(w*t)

%%%%%%-----------Modules--------------------%%%%%%%

Modules used: numpy, matplotlib, math, scipy, copy, time

%%%%%%-----------main.py--------------------%%%%%%%

main.py: This file is solely used to create a Pendulum object with various parameters
A Pendulum object has various parameters i.e if "a" is a pendulum :
a=Pendulum(mass,length,position,velocity,A,B,w,g)
All parameters for a Pendulum object are required to be an int or float type.
Some parameters if inputted as an int are automatically converted to a float type
Length and mass are required to be greater than zero also.
ValueErrors will be raised if any of these conditions are not met.

By default g is equal to 9.81 but this can be changed if needed but is not necessary to input to run.

To iterate an object "a" using the RK4 method:
a.iterateRungeKutta4(time,timestep)
"time" is the total time the object is simulated for.
"timestep" is the individual time-step for each step of the calculation. A small time-step (<1) is recommended to keep accuracy.

Time and timestep are required to be an int or float type.
Time and timestep can be negative however they must both be the same sign (else the simulation would never end).
If these conditions are not met ValueErrors are raised

The euler method is implemented into Pendulum.py however it is not very accurate and not much time was used testing it so it may not work entirely correctly.

%%%%%%-----------Example-to-run-------------%%%%%%%

a=Pendulum(1,1,math.pi/2,0,0.02,0.2,1)
a.iterateRungeKutta4(100,0.001)

%%%%%%-----------Pendulum.py----------------%%%%%%%

Pendulum.py: This file contains the calculations used to simulate the pendulum system
There should be no need to change any values inside this file.

%%%%%%----------Notes-----------------------%%%%%%%
Animation of the pendulum was planned but wasn't able to be finished in time for this version.