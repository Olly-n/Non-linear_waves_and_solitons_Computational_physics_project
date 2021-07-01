import numpy as np
# import pylab as pl
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as pl
import matplotlib.animation as animation

##### Initial condition functions

def soliton1(x, t, a):
    """basic soliton with amplitude a"""
    return 12.*(a**2.)/((np.cosh(a*((x-15.)-4*(a**2.)*t)))**2.)

def soliton2(x, t, a, b = 0.8):
    """two solitons with amplitudes a and b"""
    return 12.*(a**2.)/((np.cosh(a*((x-40.)-4*(a**2.)*t)))**2.) + 12.*(b**2.)/((np.cosh(b*((x-50.)-4*(b**2.)*t)))**2.)

def sin_1(x, t, a):
    """half period sin function"""
    if x >= 20 and x <= 70:
        return a*np.cos((x-45-t)*(2*np.pi)/100.)
    else:
        return 0

def sin_2(x, t, a):
    #sin function
    return a*np.sin(x-5/(6*np.pi))

def gauss(x, t, a):
    #gaussian shaped function
    return np.e**(-(x-25)**2)

def square(x, t, a):
    #square function
    if x>=20 and x<=30:
        return a
    else:
        return 0


##### Burgers equation Solver


class Burgers_solver:

    def __init__ (self, u, x_start, x_end, t_start, t_end, dx, dt, a, D):
        """A class designed to numerically solve the Burgers equation from an arbitrary initial condition

        params:
                u: initial conditions defined in terms of a function f(x, t, a) where a is an amplitude factor
                x_start: starting position for the x-axis
                x_end: ending position for the x-axis
                t_start: initial time coordinate
                t_end: final time coordinate
                dx: the spacial step size used during numerical integration
                dt: the temporal step size used during numerical integration
                a: amplitude of initial condition curve
                D: diffusion factor (proportional to viscosity)
        """
        self.__dx = float(dx)
        self.__dt = float(dt)
        self.__xstart = x_start
        self.__xend = x_end
        self.__tstart = t_start
        self.__tend = t_end
        self.__intfunc  = u
        self.__D = D

        h_range = np.arange(int(x_start/self.__dx), int(x_end/self.__dx)+1)*self.__dx
        self.__hrange = h_range

        t_range = np.arange(int(t_start/self.__dt),  int(t_end/self.__dt)+1)*self.__dt
        self.__trange = t_range

        points_int = []

        for i in self.__hrange:
            point = self.__intfunc(i, t_start, a)
            points_int.append(point)

        points_int = np.array(points_int)

        self.__data = [points_int]

    def plot(self, t):
        """plots the numerical solution at time t"""
        pl.figure(t)
        pl.plot(self.__hrange, self.__data[t], 'b')
        pl.xlabel("x")
        pl.ylabel("U")
        pl.show()


    def plot_one_graph(self, t_lst):
        """Plots graph of multiple solutions at times t on a single graph"""
        pl.figure(1)
        for t in t_lst:
            pl.plot(self.__hrange, self.__data[t], label="t = " + str(t*self.__dt))
        pl.legend()
        pl.xlabel("x")
        pl.ylabel("U")
        pl.show()




    def data(self):
        """Returns propagation data"""
        return self.__data


    def f(self, extra_u, i):
        """Helper function for RK4"""
        D = self.__D
        u_0 = self.__data[i]
        u_f1 = np.roll(self.__data[i]+extra_u, -1)
        u_b1 = np.roll(self.__data[i]+extra_u, 1)
        u_f2 = np.roll(self.__data[i]+extra_u, -2)
        u_b2 = np.roll(self.__data[i]+extra_u, 2)
        return (-(0.25)*(self.__dt/self.__dx)*(u_f1**2 - u_b1**2)  + D*self.__dt*(u_f1 - 2*u_0 + u_b1)/(self.__dx**2))


    def u_next(self, i):
        """RK4 implementation returning the next state"""
        u_0 = self.__data[i]
        k1 = self.f(0., i)
        k2 = self.f(0.5*k1, i)
        k3 = self.f(0.5*k2, i)
        k4 = self.f(k3, i)
        return u_0  +(1./6)*(k1 +2.*k2+2.*k3+k4)

    def propagate(self):
        """Uses the RK4 numerical method to propagate the solution through time 
            and saves the solution at each timestep.
        """
        for j in range(int(self.__dt), len(self.__trange)):
            u_j = self.u_next(j)
            self.__data.append(u_j)

    def __plot_t(self):
        """***WARNING*** FOR DEBUG ONLY. USE WITH OVER SMALL NUMBER OF TIME STEPS
        Plots a separate graph for each timestep of the solution
        """
        for i in range(len(self.__trange)):
            self.plot(i)

    def area(self):
        """Plots change in area under the solution curve with time"""
        areas = []
        for j in range(len(self.__trange)):
            strip_area = np.sum(self.__data[j])*self.__dx
            areas.append(strip_area)
        pl.figure(1)
        pl.plot(self.__trange, areas, 'b')
        pl.xlabel("t")
        pl.ylabel("Area")
        pl.show()
        return areas

    def area_2(self):
        """Plots change in area^2 under solution curve with time"""
        areas_2 = []
        for j in range(len(self.__trange)):
            strip_area = np.sum((self.__data[j]**2))*self.__dx
            areas_2.append(strip_area)
        pl.figure(1)
        pl.plot(self.__trange, areas_2, 'b')
        pl.xlabel("t")
        pl.ylabel("Area^2")
        pl.show()
        return areas_2

    def height(self, plot=True):
        """Plots change in max amplitude over time"""
        max_list = [] 
        for i in self.__data[1::]:
            max_pos = np.amax(i)
            max_list.append(max_pos)
        if plot == True:
            pl.figure(1)
            pl.plot(self.__trange, max_list)
            pl.xlabel("t")
            pl.ylabel("U_max")
            pl.show()
        return max_list

    def stab_test(self):
        """Uses height function to test the stability of the solution generated by ensuring it does not grow"""
        max_list = self.height(False)
        av_h1 = sum(max_list[0:int(len(max_list)/2)])/len(max_list[0:int(len(max_list)/2)])
        av_h2 = sum(max_list[int(len(max_list)/2):-1])/len(max_list[int(len(max_list)/2):-1])

        if av_h1*1.05 > av_h2 :
            return True
        else:
            return False

    def heatmap(self):
        """Plots heat map to show 3D dynamics of the solitons"""
        zz = np.array(self.__data)
        pl.imshow(
                np.flipud(zz.T), 
                extent=[self.__tstart, self.__tend, self.__xstart, self.__xend], 
                aspect="auto", 
                interpolation="none"
                )
        pl.ylabel("x")
        pl.xlabel("t")
        pl.colorbar(label="U")
        pl.show()

    def max(self):
        """Plots position of maximum amplitude with time"""
        max_list = [] 
        for i in self.__data[1::]:
            max_pos = np.argmax(i)
            max_list.append(max_pos)
        pl.figure(1)
        pl.plot(self.__trange, max_list)
        pl.xlabel("t")
        pl.ylabel("max_U position")
        pl.show()
        return max_list

    def vel(self):
        """Returns average velocity of the solution"""
        x_final = np.argmax(self.__data[-1])
        x_init = np.argmax(self.__data[0])
        del_x = (x_final - x_init)*self.__dx
        del_t =self.__tend-self.__tstart
        av_vel = del_x/del_t
        return av_vel

    def animation(self, show_every=1):
        """Creates an animation of the solution curve propagating.
            show_every can be used to only show every nth frame of the solution. 
            Setting show_every to 1 will animate in relative real time which can be very slow for low amplitudes.
        """
        fig = pl.figure()
        ax = pl.axes(xlim=(self.__xstart,self.__xend), ylim=(-1, np.max(self.__data[0])+1))
        pl.xlabel("x")
        pl.ylabel("U")
        line, = ax.plot([], [], lw=1)

        def init():
            return line,

        def animate(i):
            line.set_data(self.__hrange, self.__data[int(np.ceil(0.01/self.__dt))*i*show_every])
            return line,

        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(self.__data)/int(np.ceil(0.01/self.__dt)*show_every), interval=1, blit=False, repeat=True)
        pl.show()
