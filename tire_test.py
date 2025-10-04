# Import headers
from vehicle import Vehicle
from utility.parser import read_yaml
from state_models.vehicle_state import VehicleState
from scipy.optimize._root import root
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
parameters_file = 'vehicles/metric_test.yaml'
car = Vehicle(parameters_file)
vehicle_state = VehicleState(car.params)
vehicle_state.fl.fz = 800
vehicle_state.fl.gamma = np.deg2rad(5)

gamma_set = np.deg2rad(np.linspace(-15,15,5))
gamma_set2 = np.deg2rad(np.linspace(-15,15,51))
alpha_set = np.deg2rad(np.linspace(-20,20,51))
kappa_set = np.deg2rad(np.linspace(-2,2,51))

fys = np.zeros(np.size(alpha_set))
fxs = np.zeros(np.size(gamma_set2))
for gamma in gamma_set:
    vehicle_state.fl.gamma = gamma
    for i, alpha in enumerate(alpha_set):
        vehicle_state.fl.alpha = alpha
        vehicle_state.fl.eval(0.0)
        fys[i] = vehicle_state.fl.fy0
    vehicle_state.fl.alpha = 0.0
    for j, gamma in enumerate(gamma_set2):
        vehicle_state.fl.gamma = gamma
        vehicle_state.fl.eval(1.0)
        fxs[j] = vehicle_state.fl.fx0
    # plt.plot(alpha_set,fys)
    print(len(fxs))

    plt.plot(gamma_set2,fxs)

plt.show()



