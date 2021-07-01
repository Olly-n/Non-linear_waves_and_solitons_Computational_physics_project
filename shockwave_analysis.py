from shockwave_solver import *

""" Use this scrip to model the burges equation"""
# init_curve = soliton1
# init_curve = soliton2
# init_curve = sin_1
# init_curve = sin_2
# init_curve = gauss
init_curve = square

z = Burgers_solver(init_curve, 0, 50., 0,12, 0.1, 0.001, 2.5, 0.1)
z.propagate()

z.animation(10)

z.plot_one_graph([0,3000,8000,12000])
z.area()
z.area_2()
# z.height()
z.heatmap()