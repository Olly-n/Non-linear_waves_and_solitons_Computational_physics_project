from KDV_solver_and_stability import *

"""Use this script to run different aspects of the KDV analysis solver and plot various graphs"""

##### KDV solver and animations

# init_curve = soliton1
# init_curve = soliton2
# init_curve = sin
# init_curve = gaussian
init_curve = gauss_plus_soliton
# init_curve = square

soln = KDV_solver(init_curve, 0, 100., 0, 50., 0.1, 0.001, 3)
soln.propagate()

soln.animation(20)
soln.plot_one_graph([0,3000,10000])
soln.heatmap()


##### KDV solver analysis

# velocity(0.1,7,0.1)

# find_h_crit(1,10,0.5,0.08,1.,0.005)
