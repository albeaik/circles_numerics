classdef UserDefinedCoupling < handle
    %USERDEFINEDCOUPLING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        activate_coupling;
        
        Autonomous_state;
        Autonomous_control_u;           %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now
        time;
        time_step_number;
        
        Autonomous_state_history;
    end
    
    methods
        function obj = UserDefinedCoupling(activate_coupling, initial_autonomous_state, autonomous_control_u)
            %default constructor object
            
            obj = obj@handle();
            
            obj.activate_coupling = activate_coupling;        
            
            obj.Autonomous_state = initial_autonomous_state;  %initial state of the autonomous cars.. each row is [position, speed]
            obj.Autonomous_control_u = autonomous_control_u;           %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now

            obj.time = 0;
            obj.time_step_number = 1;
            
            obj.Autonomous_state_history{1} = obj.Autonomous_state;
        end
        
        function SimulateStep(obj, dt, PDEModel, PDE_evolving_mesh_obj)
            if(obj.activate_coupling)
                solution_mesh = PDE_evolving_mesh_obj;

                Autonomous_state_t = obj.Autonomous_state;
                Autonomous_control_u = obj.Autonomous_control_u(:, mod(obj.time_step_number - 1, size(obj.Autonomous_control_u, 2)) + 1);

                % simulate the autonomous car ODE
                Autonomous_state_tau = Autonomous_state_t + dt*[PDEModel.G(solution_mesh, obj, Autonomous_state_t), ...
                                                PDEModel.H(solution_mesh, obj, Autonomous_state_t) + Autonomous_control_u];

                % update object
                obj.time = obj.time + dt;
                obj.time_step_number = obj.time_step_number + 1;
                obj.Autonomous_state = Autonomous_state_tau;
                obj.Autonomous_state_history{obj.time_step_number} = obj.Autonomous_state;
            end
        end
        
        function VisualizeStep(obj, tstep)
            if(obj.activate_coupling)
                if(tstep < 0) %visualize last time step if not given or negative values
                    tstep = obj.time_step_number;
                end

                Autonomous_state_history = obj.Autonomous_state_history;
                hold on
                plot3(Autonomous_state_history{tstep}(:, 1), Autonomous_state_history{tstep}(:, 2), 10*ones(size(Autonomous_state_history{tstep}(:, 2))), 'r*')
                hold off
            end
        end
    end
    
end

