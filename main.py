import math
from Pendulum import Pendulum

#a=Pendulum(10,10,math.pi/2,0,0.02,0.2,1)    ##Forced damped non-linear pendulum
#a.iterateRungeKutta4(100,0.001)

#b=Pendulum(10,10,math.pi/2,0,0,0.02,0)      ##Free Pendulum
#b.iterateRungeKutta4(100,0.001)

c=Pendulum(10,10,math.pi/2,0,0,0,0)         ##simple pendulum
c.iterateRungeKutta4(100,0.001)

print("Done!")