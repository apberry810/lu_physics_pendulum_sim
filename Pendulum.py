import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import copy
import time as mytime

class Pendulum:
    """
    Class to model a forced pendulum system with different parameters.
    Default values of acceleration due to gravity, length,
    The equation of the pendulum is of the form:
    d^2(theta)/dt^2 + B*d(theta)/dt+g*sin(theta)/l=A*cos(w*t)
    Units:
    Angles - Radians (rad)
    Masses - Kilograms (kg)
    Lengths - Metres (m)
    Accelerations - radians per second per second (rads^-2)
    Velocities - radians per second (rads^-1)
    Positions - radians from vertical (rad)
    Energies - Joules (J)
    Momenta - kilogram metre per second (kgms^-1)
    Angular Momenta - kilogram metres^2 per second (kgm^2 s^-1)
    """
    #G=9.81 #ms^-2
    start_time=mytime.time()
    mass=1.0
    position=math.pi
    velocity=0.0
    A=0.0
    B=0.0
    w=0.0
    length=1.0
    g=9.81

    def __init__(self,mass,length,position,velocity,A,B,w,g=9.81):
        if not isinstance(mass,(int,float)):
            raise ValueError("mass is not an int or a float!")
        if not mass>0:
            raise ValueError("mass is not greater than 0!")        
        if not isinstance(length,(int,float)):
            raise ValueError("length is not an int or a float!")
        if not length>=0:
            raise ValueError("length is not greater or equal to 0!")        
        if not isinstance(position,(int,float)):
            raise ValueError("position is not an int or float!")
        if isinstance(position,int):
            position=float(position)
        if not isinstance(velocity,(int,float)):
            raise ValueError("velocity is not an int or float!")        
        if isinstance(velocity,int):
            velocity=float(velocity)
        if not isinstance(A,(int,float)):
            raise ValueError("A is not an int or a float!")
        if isinstance(A, int):
            A=float(A)
        if not isinstance(B,(int,float)):
            raise ValueError("B is not an int or a float!")
        if isinstance(B, int):
            B=float(B)
        if not isinstance(w,(int,float)):
            raise ValueError("w is not an int or a float!")
        if isinstance(w, int):
            w=float(w)
        if not isinstance(g,(int,float)):
            raise ValueError("g is not an int or a float!")
        if isinstance(g,int):
            g=float(g)
        self.mass=mass
        self.length=length
        self.position=position
        self.velocity=velocity
        self.A=A
        self.B=B
        self.w=w
        self.g=9.81
        super().__init__()

    def __repr__(self):
        return 'Mass: %f, Length: %f, Position: %.6f, Velocity %.6f, A: %.6f, B: %.6f, w: %.6f, g: %f'%(self.mass,self.length,self.position,self.velocity,self.A,self.B,self.w,self.g)#super().__repr__()

    @staticmethod
    def printcoeffs(self):
        return "A=", self.A, "B=", self.B, "w=", self.w, "g=", self.g, "m=", self.mass, "length=", self.length, "position(rads)=", self.position, "velocity(rads^-1)=", self.velocity

    @staticmethod
    def calculateMomentum(self):
        """Calculates momentum vector of solid mass in cartesian coordinates (x,y)"""
        momentumarray=np.zeros((2))
        momentumarray[0]=self.length*math.cos(self.velocity)*self.mass
        momentumarray[1]=self.length*math.sin(self.velocity)*self.mass
        return momentumarray

    @staticmethod
    def calculateAngularMomentum(self):
        """Calculates angular momentum of solid mass"""
        return self.velocity*self.mass*self.length

    @staticmethod
    def calculateKineticEnergy(self):
        """Calculates kinetic energy of the system (solid mass)"""       
        return 0.5*self.mass*self.velocity**2

    @staticmethod
    def calculateGPE(self):
        """Calculates current gravitational potential energy"""
        return self.mass*self.g*self.length*(1-math.cos(self.position))

    @staticmethod
    def calculateTotalEnergy(self):
        """Calculates total energy of the system (solid mass)"""
        kineticenergy=self.calculateKineticEnergy(self)
        potentialenergy=self.calculateGPE(self)
        return kineticenergy+potentialenergy

    @staticmethod
    def calculateLagrangian(self):
        """Calculates Lagrangian of the system (solid mass)"""
        kineticenergy=self.calculateKineticEnergy(self)
        potentialenerg=self.calculateGPE(self)
        return kineticenergy-potentialenerg

    @staticmethod
    def getcartesiancoords(self):
        """Calculates 2D cartesian coordinates in an array"""
        cartesian=np.zeros((2))
        cartesian[0]=self.length*math.cos(self.position)
        cartesian[1]=self.length*math.sin(self.position)
        return cartesian

    @staticmethod
    def getcartesianvelocity(self):
        """Calcaulates 2D cartesian velocity vector in an array"""
        cartesianvelocity=np.zeros((2))
        cartesianvelocity[0]=self.length*math.cos(self.velocity)
        cartesianvelocity[1]=self.length*math.sin(self.velocity)
        return cartesianvelocity

    def iterateEuler(self, time, timestep, type):
        """
        Using a second-order modification of the basic Euler numerical method iterates a Pendulum object for a given time and timestep.
        time: Total length of time the simulation is simulated for.
        timestep: size of the the timestep used in the algorithms. A smaller timestep results in larger precision but takes longer.
        """
        totalT=0.0
        omega_n=self.velocity
        theta_n=self.position
        g=self.g
        A=self.A
        B=self.B
        l=self.length
        w=self.w
        ##arrays created for data collection
        positionarray=np.array([])          ##1
        velocityarray=np.array([])          ##2
        data=np.array([])                   ##3, time array
        totalenergyarray=np.array([])       ##4/entered
        kineticenergyarray=np.array([])     ##5/entered
        potentialenergyarray=np.array([])   ##6/entered
        lagrangianarray=np.array([])        ##7/entered
        momentumarray=np.array([])          ##8/entered
        angularmomentumarray=np.array([])   ##9/
        cartesianpositionarray=np.array([]) ##10    
        while(totalT<time):
            omega_n+=timestep*(-self.B*omega_n-self.g*math.sin(theta_n)/self.length+self.A*math.cos(self.w*totalT))
            theta_n+=timestep*omega_n
            self.velocity=omega_n
            self.position=theta_n
            #record theta then loop above should be correct
            totalT+=timestep
        print("Final time:",max(data))
        print("This RK4 simulation took", mytime.time()-self.start_time, "seconds from start to exit")
        fig=plt.figure()
        #position and velocity against time
        plt.subplot()
        pos=fig.gca()
        pos.plot(data,positionarray,color="r",linestyle="solid",label="position (polar coordinates)",linewidth=1.5)
        pos.plot(data,velocityarray,color="b",linestyle="solid",label="velocity (polar coordinates)",linewidth=1.5)
        pos.set_title("Position and Velocity against Time")
        pos.set_xlabel("time (sec)")
        pos.set_ylabel("position (rad)")
        pos.legend()
        pos.grid()
        plt.show()
        #energies against time
        plt.subplot()
        plt.plot(data,totalenergyarray,color="r",linestyle=(0,(3,1,1,1)),label="Total energy",linewidth=1.5)
        plt.plot(data,lagrangianarray,color="b",linestyle=(0,(5,1)),label="Lagrangian of system",linewidth=1.5)
        plt.plot(data,kineticenergyarray,color="g",linestyle="solid",label="Kinetic energy",linewidth=1.5)
        plt.plot(data,potentialenergyarray,color="m",linestyle="solid",label="Gravitational potential energy",linewidth=1.5)
        plt.title("Energies against Time")
        plt.xlabel("time (s)")
        plt.ylabel("Energy (J)")
        plt.legend()
        plt.grid()
        plt.show()
        #momenta against time
        plt.subplot()
        plt.plot(data,momentumarray[:,0],color="r",linestyle="solid",label="Momentum in x-direction",linewidth=1.5)
        plt.plot(data,momentumarray[:,1],color="b",linestyle="solid",label="Momentum in y-direction",linewidth=1.5)
        plt.title("Momentum in x and y Directions against Time")
        plt.xlabel("time (s)")
        plt.ylabel("Momentum (kgm/s)")
        plt.legend()
        plt.grid()
        plt.show()
        #angmomenta against time
        plt.subplot()
        plt.plot(data,angularmomentumarray,color="r",linestyle="solid",label="Angular momentum",linewidth=1.5)
        plt.title("Angular Momentum against Time")
        plt.xlabel("time (s)")
        plt.ylabel("Angular momentum (kgm^2/s)")
        plt.legend()
        plt.grid()
        plt.show()
        return    

    def iterateRungeKutta4(self, time, timestep):
        """
        Using a second order modification of the Runge-Kutta 4 (RK4) method iterates a Pendulum object for a given timestep for a given length of time.
        time: Total length of time the simulation is simulated for.
        timestep: size of the the timestep used in the algorithms. A smaller timestep results in larger precision but takes longer.
        """
        if not isinstance(time,(int,float)):
            raise ValueError("time is not an int or float!")
        if not isinstance(timestep,(int,float)):
            raise ValueError("timestep is not an int or float!")
        if (timestep==0):
            raise ValueError("timestep cannot be zero!")
        if (time*timestep<0):
            raise ValueError("Time or timestep is negative while the other is positive. Both values must have the same sign!")
        print("RK4 iteration started for",time,"seconds with a timestep of",timestep,"seconds for a total of",time/timestep,"steps")
        totalT=0.0
        omega_n=self.velocity
        theta_n=self.position
        g=self.g
        A=self.A
        B=self.B
        l=self.length
        w=self.w
        ##arrays created for data collection
        positionarray=np.array([])          ##1/entered
        velocityarray=np.array([])          ##2/entered
        data=np.array([])                   ##3
        totalenergyarray=np.array([])       ##4/entered
        kineticenergyarray=np.array([])     ##5/entered
        potentialenergyarray=np.array([])   ##6/entered
        lagrangianarray=np.array([])        ##7/entered
        momentumarray=np.array([])          ##8/entered
        angularmomentumarray=np.array([])   ##9/entered
        cartesianpositionarray=np.array([]) ##10/entered
        while(totalT<time):
            k1theta=omega_n
            k1omega=-B*omega_n-g*math.sin(theta_n)/l+A*math.cos(w*totalT)
            k2theta=omega_n+timestep*k1omega/2
            k2omega=-B*(omega_n+timestep*k1omega/2)-g*math.sin(theta_n+timestep*k1omega/2)/l+A*math.cos(w*(totalT+timestep/2))
            k3theta=omega_n+timestep*k2omega/2
            k3omega=-B*(omega_n+timestep*k2omega/2)-g*math.sin(theta_n+timestep*k2omega/2)/l+A*math.cos(w*(totalT+timestep/2))
            k4theta=omega_n+timestep*k3omega
            k4omega=-B*(omega_n+timestep*k3omega)-g*math.sin(theta_n+timestep*k3theta)/l+A*math.cos(w*(totalT+timestep))
            theta_n+=timestep*(k1theta+2*k2theta+2*k3theta+k4theta)/6
            omega_n+=timestep*(k1omega+2*k2omega+2*k3omega+k4omega)/6            
            self.velocity=omega_n
            self.position=theta_n
            totalT+=timestep
            ####arrays update below
            velocityarray=np.append(velocityarray,omega_n)      
            positionarray=np.append(positionarray,theta_n)
            data=np.append(data,totalT)
            totalenergyarray=np.append(totalenergyarray,self.calculateTotalEnergy(self))
            kineticenergyarray=np.append(kineticenergyarray,self.calculateKineticEnergy(self))
            potentialenergyarray=np.append(potentialenergyarray,self.calculateGPE(self))
            lagrangianarray=np.append(lagrangianarray,self.calculateLagrangian(self))
            if not momentumarray.size==0:
                momentumarray=np.vstack((momentumarray,self.calculateMomentum(self)))
            if momentumarray.size==0:
                momentumarray=np.append(momentumarray,self.calculateMomentum(self))
            angularmomentumarray=np.append(angularmomentumarray,self.calculateAngularMomentum(self))
            newcartesianposition=np.array([self.length*math.cos(theta_n),self.length*math.sin(theta_n)])
            if not cartesianpositionarray.size==0:
                cartesianpositionarray=np.vstack((cartesianpositionarray,newcartesianposition))
            if cartesianpositionarray.size==0:
                cartesianpositionarray=np.append(cartesianpositionarray,newcartesianposition)
        print("Final time:",max(data),"seconds")
        print("This RK4 simulation took", mytime.time()-self.start_time, "seconds from start to exit")
        fig=plt.figure()
        #graph - position and velocity against time
        plt.subplot()
        pos=fig.gca()
        pos.plot(data,positionarray,color="r",linestyle="solid",label="position (polar coordinates)",linewidth=1.5)
        pos.plot(data,velocityarray,color="b",linestyle="solid",label="velocity (polar coordinates)",linewidth=1.5)
        pos.set_title("Position and Velocity against Time")
        pos.set_xlabel("time (sec)")
        pos.set_ylabel("position (rad)/velocity (rad/s)")
        pos.legend()
        pos.grid()
        plt.show()
        #graph - energies against time
        plt.subplot()
        plt.plot(data,totalenergyarray,color="r",linestyle=(0,(3,1,1,1)),label="Total energy",linewidth=1.5)
        plt.plot(data,lagrangianarray,color="b",linestyle=(0,(5,1)),label="Lagrangian of system",linewidth=1.5)
        plt.plot(data,kineticenergyarray,color="g",linestyle="solid",label="Kinetic energy",linewidth=1.5)
        plt.plot(data,potentialenergyarray,color="m",linestyle="solid",label="Gravitational potential energy",linewidth=1.5)
        plt.title("Energies against Time")
        plt.xlabel("time (s)")
        plt.ylabel("Energy (J)")
        plt.legend()
        plt.grid()
        plt.show()
        #graph - momenta against time
        plt.subplot()
        plt.plot(data,momentumarray[:,0],color="r",linestyle="solid",label="Momentum in x-direction",linewidth=1.5)
        plt.plot(data,momentumarray[:,1],color="b",linestyle="solid",label="Momentum in y-direction",linewidth=1.5)
        plt.title("Momentum in x and y Directions against Time")
        plt.xlabel("time (s)")
        plt.ylabel("Momentum (kgm/s)")
        plt.legend()
        plt.grid()
        plt.show()
        #graph - angmomenta against time
        plt.subplot()
        plt.plot(data,angularmomentumarray,color="r",linestyle="solid",label="Angular momentum",linewidth=1.5)
        plt.title("Angular Momentum against Time")
        plt.xlabel("time (s)")
        plt.ylabel("Angular momentum (kgm^2/s)")
        plt.legend()
        plt.grid()
        plt.show()
        return