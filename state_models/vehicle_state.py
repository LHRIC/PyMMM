from state_models.dynamic_module import DynModel
from state_models.powertrain_module import PtnModel
from state_models.aero_module import AeroModel
from state_models.tire_state import TireState
from scipy.constants import g
import numpy as np

class VehicleState:
    def __init__(self,params:dict) -> None:
        # Pull parameter dictionary
        self.params = params['params']
        # print(self.params)
        # Initialize models
        self.dyn = DynModel(self)
        self.ptn = PtnModel(self)
        self.aero = AeroModel(self)
        # Initialize Tire models
        self.fl = TireState(self)
        self.fr = TireState(self)
        self.rl = TireState(self)
        self.rr = TireState(self)
        # Sprung mass state
        self.roll = 0
        self.pitch = 0



    def _residuals(self):
        mz = np.linalg.matmul(self.params['inertia_tensor'],[0,0,self.psi_ddt])
        self.moments.append(mz)
        # print (f'Forces: {self.forces}')
        # print (f'Moments: {self.moments}')
        sum_forces = np.add.reduce(self.forces)
        sum_moments = np.add.reduce(self.moments)
        # print (f'sum_Forces: {sum_forces}')
        # print (f'sum_Moments: {sum_moments}')
        residuals = [sum_forces[0], sum_forces[1], sum_moments[2]]
        return residuals

    def eval(self,v,psi_dt,beta,delta,eta,x_ddt,y_ddt,psi_ddt,residuals):
        # Initialize Tires
        for tire in ['fl','fr','rl','rr']:
            getattr(self,tire).fz = 0
        self.forces=[]
        self.moments=[]
        self.v=v # Tangential velocity
        self.psi_dt = psi_dt # Yaw rate
        self.beta=beta # Body slip
        self.delta=delta # Steered angle
        self.eta=eta # 'Throttle' 
        self.x_ddt=x_ddt # Ax
        self.y_ddt=y_ddt # Ay
        # print(f'Ax {self.x_ddt} Ay {self.y_ddt}')
        self.psi_ddt=psi_ddt # Yaw accel
        # Each function updates: self.[fl,fr,rl,rr] | self.forces | self.moments
        self.dyn.static_weight(self) # Gravity effects
        self.dyn.weight_transfer(self) # Acceleration effects
        self.aero.downforce(self) # Downforce
        self.dyn.kinematic_eval(self) # Suspension travel and inclination angle
        self.ptn.torque(self) #Powertrain and braking, updates max Fx for tires
        self.dyn.steering(self) # Evaluate steering angles
        self.fl.eval(self.eta) # Evaluate tire forces
        self.fr.eval(self.eta)
        self.rl.eval(self.eta)
        self.rr.eval(self.eta)
        self.dyn.tire_forces(self) # Apply tire forces to car
        # if np.isnan(self.fl.fz) or np.isnan(self.fr.fz):
            # print(f'NAN @ Ax {self.x_ddt}, Ay {self.y_ddt}, Psi_dt {self.psi_dt}')
        if residuals == True:
            r = self._residuals() 
            if any(np.isnan(r)):
                print(f'NAN @')
                print(f'FL: {self.fl.f_vec}')
                print(f'FR: {self.fr.f_vec}')
                print(f'RL: {self.rl.f_vec}')
                print(f'RR: {self.rr.f_vec}')
            return r
        

