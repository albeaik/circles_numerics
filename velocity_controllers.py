"""Contains a list of custom velocity controllers."""

from flow.controllers.base_controller import BaseController
import numpy as np


class FollowerStopper(BaseController):
    """Inspired by Dan Work's... work.

    Dissipation of stop-and-go waves via control of autonomous vehicles:
    Field experiments https://arxiv.org/abs/1705.01693

    Usage
    -----
    See base class for example.

    Parameters
    ----------
    veh_id : str
        unique vehicle identifier
    v_des : float, optional
        desired speed of the vehicles (m/s)
    """

    def __init__(self,
                 veh_id,
                 car_following_params,
                 v_des=15,
                 danger_edges=None,
                 control_length=None,
                 no_control_edges=None):
        """Instantiate FollowerStopper."""
        BaseController.__init__(
            self, veh_id, car_following_params, delay=0.0,
            fail_safe='safe_velocity')

        # desired speed of the vehicle
        self.v_des = v_des

        # maximum achievable acceleration by the vehicle
        self.max_accel = car_following_params.controller_params['accel']

        # other parameters
        self.dx_1_0 = 4.5
        self.dx_2_0 = 5.25
        self.dx_3_0 = 6.0
        self.d_1 = 1.5
        self.d_2 = 1.0
        self.d_3 = 0.5

        self.danger_edges = danger_edges if danger_edges else {}
        self.control_length = control_length
        self.no_control_edges = no_control_edges

    def find_intersection_dist(self, env):
        """Find distance to intersection.

        Parameters
        ----------
        env : flow.envs.Env
            see flow/envs/base.py

        Returns
        -------
        float
            distance from the vehicle's current position to the position of the
            node it is heading toward.
        """
        edge_id = env.k.vehicle.get_edge(self.veh_id)
        # FIXME this might not be the best way of handling this
        if edge_id == "":
            return -10
        if 'center' in edge_id:
            return 0
        edge_len = env.k.network.edge_length(edge_id)
        relative_pos = env.k.vehicle.get_position(self.veh_id)
        dist = edge_len - relative_pos
        return dist

    def get_accel(self, env):
        """See parent class."""
        if env.time_counter < env.env_params.warmup_steps * env.env_params.sims_per_step:
            return None
        else:
            lead_id = env.k.vehicle.get_leader(self.veh_id)
            this_vel = env.k.vehicle.get_speed(self.veh_id)
            lead_vel = env.k.vehicle.get_speed(lead_id)

            if self.v_des is None:
                return None

            if lead_id is None:
                v_cmd = self.v_des
            else:
                dx = env.k.vehicle.get_headway(self.veh_id)
                dv_minus = min(lead_vel - this_vel, 0)

                dx_1 = self.dx_1_0 + 1 / (2 * self.d_1) * dv_minus**2
                dx_2 = self.dx_2_0 + 1 / (2 * self.d_2) * dv_minus**2
                dx_3 = self.dx_3_0 + 1 / (2 * self.d_3) * dv_minus**2
                v = min(max(lead_vel, 0), self.v_des)
                # compute the desired velocity
                if dx <= dx_1:
                    v_cmd = 0
                elif dx <= dx_2:
                    v_cmd = v * (dx - dx_1) / (dx_2 - dx_1)
                elif dx <= dx_3:
                    v_cmd = v + (self.v_des - this_vel) * (dx - dx_2) \
                            / (dx_3 - dx_2)
                else:
                    v_cmd = self.v_des

            edge = env.k.vehicle.get_edge(self.veh_id)

            if edge == "":
                return None

            if (self.find_intersection_dist(env) <= 10 and
                    env.k.vehicle.get_edge(self.veh_id) in self.danger_edges) or \
                    env.k.vehicle.get_edge(self.veh_id)[0] == ":" \
                    or (self.control_length and (env.k.vehicle.get_x_by_id(self.veh_id) < self.control_length[0]
                        or env.k.vehicle.get_x_by_id(self.veh_id) > self.control_length[1])) \
                    or (self.no_control_edges is not None and len(self.no_control_edges) > 0
                        and edge in self.no_control_edges):
                return None
            else:
                # compute the acceleration from the desired velocity
                return np.clip((v_cmd - this_vel) / env.sim_step, -np.abs(self.max_deaccel), self.max_accel)


class NonLocalFollowerStopper(FollowerStopper):
    """Follower stopper that uses the average system speed to compute its acceleration."""

    def get_accel(self, env):
        """See parent class."""
        lead_id = env.k.vehicle.get_leader(self.veh_id)
        this_vel = env.k.vehicle.get_speed(self.veh_id)
        lead_vel = env.k.vehicle.get_speed(lead_id)
        self.v_des = np.mean(env.k.vehicle.get_speed(env.k.vehicle.get_ids()))

        if self.v_des is None:
            return None

        if lead_id is None:
            v_cmd = self.v_des
        else:
            dx = env.k.vehicle.get_headway(self.veh_id)
            dv_minus = min(lead_vel - this_vel, 0)

            dx_1 = self.dx_1_0 + 1 / (2 * self.d_1) * dv_minus ** 2
            dx_2 = self.dx_2_0 + 1 / (2 * self.d_2) * dv_minus ** 2
            dx_3 = self.dx_3_0 + 1 / (2 * self.d_3) * dv_minus ** 2
            v = min(max(lead_vel, 0), self.v_des)
            # compute the desired velocity
            if dx <= dx_1:
                v_cmd = 0
            elif dx <= dx_2:
                v_cmd = v * (dx - dx_1) / (dx_2 - dx_1)
            elif dx <= dx_3:
                v_cmd = v + (self.v_des - this_vel) * (dx - dx_2) \
                        / (dx_3 - dx_2)
            else:
                v_cmd = self.v_des

        edge = env.k.vehicle.get_edge(self.veh_id)

        if edge == "":
            return None
        else:
            # compute the acceleration from the desired velocity
            return (v_cmd - this_vel) / env.sim_step


class PISaturation(BaseController):
    """Inspired by Dan Work's... work.

    Dissipation of stop-and-go waves via control of autonomous vehicles:
    Field experiments https://arxiv.org/abs/1705.01693

    Usage
    -----
    See base class for example.

    Parameters
    ----------
    veh_id : str
        unique vehicle identifier
    car_following_params : flow.core.params.SumoCarFollowingParams
        object defining sumo-specific car-following parameters
    """

    def __init__(self, veh_id, car_following_params):
        """Instantiate PISaturation."""
        BaseController.__init__(self, veh_id, car_following_params, delay=0.0)

        # maximum achievable acceleration by the vehicle
        self.max_accel = car_following_params.controller_params['accel']

        # history used to determine AV desired velocity
        self.v_history = []

        # other parameters
        self.gamma = 2
        self.g_l = 7
        self.g_u = 30
        self.v_catch = 1

        # values that are updated by using their old information
        self.alpha = 0
        self.beta = 1 - 0.5 * self.alpha
        self.U = 0
        self.v_target = 0
        self.v_cmd = 0

    def get_accel(self, env):
        """See parent class."""
        lead_id = env.k.vehicle.get_leader(self.veh_id)
        lead_vel = env.k.vehicle.get_speed(lead_id)
        this_vel = env.k.vehicle.get_speed(self.veh_id)

        dx = env.k.vehicle.get_headway(self.veh_id)
        dv = lead_vel - this_vel
        dx_s = max(2 * dv, 4)

        # update the AV's velocity history
        self.v_history.append(this_vel)

        if len(self.v_history) == int(38 / env.sim_step):
            del self.v_history[0]

        # update desired velocity values
        v_des = np.mean(self.v_history)
        v_target = v_des + self.v_catch \
            * min(max((dx - self.g_l) / (self.g_u - self.g_l), 0), 1)

        # update the alpha and beta values
        alpha = min(max((dx - dx_s) / self.gamma, 0), 1)
        beta = 1 - 0.5 * alpha

        # compute desired velocity
        self.v_cmd = beta * (alpha * v_target + (1 - alpha) * lead_vel) \
            + (1 - beta) * self.v_cmd

        # compute the acceleration
        accel = (self.v_cmd - this_vel) / env.sim_step

        return min(accel, self.max_accel)
    
    
class MicromodelController(BaseController):

    def __init__(self,
                 veh_id,
                 car_following_params,
                 v_des=13.336,
                 Penetration_rate=0.05,
                 danger_edges=None,
                 control_length=None,
                 no_control_edges=None):
        """Instantiate FollowerStopper."""
        BaseController.__init__(
            self, veh_id, car_following_params, delay=1.0,
            fail_safe='safe_velocity')
    
        # desired speed of the vehicle
    #    self.v_des = v_des
    #    self.entry_headway = solve TRAFFIC_FLOW
    #    self.TRAFFIC_FLOW = solve TRAFFIC_FLOW
    #Design parameters
        self.d0 = 2.5
        self.a = 0.2
        self.b = 0.001
        self.v_max = 30
        #self.TRAFFIC_INFLOW_speed = 11
 #       if env == None:
        self.v_des = v_des
        self.v_fast = 24
 #       else:            
  #          self.v_des = min(self.a+self.b*env.time_counter * env.sim_step,1)*v_des
#        self.v_des = min(self.a+self.b*env.time_counter * env.sim_step,1)*TRAFFIC_INFLOW_speed
#        self.v_des = min(self.a+self.b*env.time_counter * env.sim_step,1)*v_max*((np.tanh(env.entry_headway/env.h_st-2)+np.tanh(2))/(1+np.tanh(2)))
        self.k = 1.0
        self.dx_1_0 = 4.5
        self.dx_2_0 = 5.25
        self.dx_3_0 = 6.0
        self.d_1 = 1.5
        self.d_2 = 1.0
        self.d_3 = 0.5
        self.d = 1
        self.ref_time = 0
        # maximum achievable acceleration by the vehicle
        self.max_accel = car_following_params.controller_params['accel']
    
    
        # other parameters
        self.danger_edges = danger_edges if danger_edges else {}
    
    #def ov_model(h,self,env):
    #    return self.v_max*((np.tanh(h/env.h_st-2)+np.tanh(2))/(1+np.tanh(2)))/h-3600/self.TRAFFIC_FLOW

   # def ov_model(h):
   #     return v_max*((np.tanh(h/h_st-2)+np.tanh(2))/(1+np.tanh(2)))*3600/(h*TRAFFIC_FLOW)-1

    
#    self.entry_hd = scipy.optimize.fsolve(ov_model,11,(self,env))    
    
    def find_intersection_dist(self, env):
        """Find distance to intersection.
    
        Parameters
        ----------
        env : flow.envs.Env
            see flow/envs/base.py
    
        Returns
        -------
        float
            distance from the vehicle's current position to the position of the
            node it is heading toward.
        """
        edge_id = env.k.vehicle.get_edge(self.veh_id)
        # FIXME this might not be the best way of handling this
        if edge_id == "":
            return -10
        if 'center' in edge_id:
            return 0
        edge_len = env.k.network.edge_length(edge_id)
        relative_pos = env.k.vehicle.get_position(self.veh_id)
        dist = edge_len - relative_pos
        return dist
    
    def get_accel(self, env):
        """See parent class."""
        lead_id = env.k.vehicle.get_leader(self.veh_id)
        this_vel = env.k.vehicle.get_speed(self.veh_id)
        lead_vel = env.k.vehicle.get_speed(lead_id)
        desired_lane = env.k.vehicle.get_lane(self.veh_id)
    
        veh_ids = env.k.vehicle.get_ids()
#        lanes = env.k.vehicle.get_lanes()
        in_lane = veh_ids
        #[veh_id for veh_id in veh_ids if env.k.vehicle.get_lane(veh_id) == desired_lane]
        avg_speed = np.nan_to_num(np.mean(env.k.vehicle.get_speed(in_lane)))

#        s = 0
#        for _veh_id_v in ...:
#            s = s + env.k.vehicle.get_speed(_veh_id)
#        while d ==0:
#            s = s + env.k.vehicle.get_speed(_veh_id)
        
#        avg_speed = speeds(this_lane)
        if env.time_counter < 1500:
            v_des2 = self.v_fast
            
            #v_des2 = avg_speed
        if self.d == 0:
#        elif self.d == 0:
            v_des2 = self.v_des*min(0.8+0.01*(env.time_counter * env.sim_step - self.ref_time),1)
        else:
            v_des2 = self.v_des 
            #v_des2 = self.v_des * min(max(avg_speed,self.a + self.b * env.time_counter * env.sim_step), 1)
        if v_des2 is None:
            return None
    
        if lead_id is None:
            v_cmd = self.v_des
            
        else:
            dx = env.k.vehicle.get_headway(self.veh_id)
            dv_minus = min(lead_vel - this_vel, 0)
            
            dx_1 = self.dx_1_0 + 1 / (2 * self.d_1) * dv_minus**2
            dx_2 = self.dx_2_0 + 1 / (2 * self.d_2) * dv_minus**2
            dx_3 = self.dx_3_0 + 1 / (2 * self.d_3) * dv_minus**2
            
            v = min(max(lead_vel, 0), v_des2)
            # compute the desired velocity
            
#            if dx <= dx_1:
#                v_cmd = 0
#            elif dx <= dx_2:
#                v_cmd = v * (dx - dx_1) / (dx_2 - dx_1)
#            elif dx <= dx_3:
#                v_cmd = v + (self.v_des - this_vel) * (dx - dx_2) \
#                        / (dx_3 - dx_2)
#            else:
            if dx >= self.dx_3_0*7:
                v_cmd = v_des2
#                v_cmd = self.v_fast*min(1,(dx-self.dx_3_0*7))+v_des2*(1-(dx-self.dx_3_0*7))
                self.d = 1
            elif dx <= self.dx_1_0:
                v_cmd = v
                self.d = 0
                self.ref_time = env.time_counter * env.sim_step
            else:
                v_cmd = v_des2
            
            if env.k.vehicle.get_distance(self.veh_id) < 100:
                v_cmd = self.v_fast
    
#Faire en sorte que quand on ralenti on change de setting pour la vitesse pour la rendre plus faible
#Ajuster en fonction de la zone pour les differentes simulations
#Trouver un moyen de garder le passé en mémoire pour avoir accès à un avg speed cohérent.                
        edge = env.k.vehicle.get_edge(self.veh_id)
    
        if edge == "":
            return None
    
        if self.find_intersection_dist(env) <= 10 and \
                env.k.vehicle.get_edge(self.veh_id) in self.danger_edges or \
                env.k.vehicle.get_edge(self.veh_id)[0] == ":":
            return None
        else:
            # compute the acceleration from the desired velocity
            return self.k*(v_cmd - this_vel) / env.sim_step
            # return k*(v_cmd - this_vel)
