classdef UserDefinedCoupling < handle
    %USERDEFINEDCOUPLING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Autonomous_state;
        %Autonomous_init_state;  %initial state of the autonomous cars.. each row is [position, speed]
        Autonomous_control_u;           %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now
        time;
        time_step_number;
        
        Autonomous_state_history;
    end
    
    methods
        function InitializeDefault(obj)
            obj.Autonomous_state = [7, 7; 9, 7];  %initial state of the autonomous cars.. each row is [position, speed]
            obj.Autonomous_control_u = [-2; -2];           %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now

            obj.time = 0;
            obj.time_step_number = 1;
            
            obj.Autonomous_state_history{1} = obj.Autonomous_state;
        end
        
        function SimulateStep(obj, dt, PDE_evolving_mesh_obj)
            solution_mesh = PDE_evolving_mesh_obj;
            
            Autonomous_state_t = obj.Autonomous_state;
            Autonomous_control_u = obj.Autonomous_control_u;
            
            % simulate the autonomous car ODE
            Autonomous_state_tau = Autonomous_state_t + dt*[Autonomous_state_t(:, 2), (H_IDM(solution_mesh, obj, Autonomous_state_t) + Autonomous_control_u)];
            
            % update object
            obj.time = obj.time + dt;
            obj.time_step_number = obj.time_step_number + 1;
            obj.Autonomous_state = Autonomous_state_tau;
            obj.Autonomous_state_history{obj.time_step_number} = obj.Autonomous_state;
        end
        
        function VisualizeStep(obj)
            Autonomous_state_history = obj.Autonomous_state_history;
            hold on
            plot3(Autonomous_state_history{end}(:, 1), Autonomous_state_history{end}(:, 2), 10*ones(size(Autonomous_state_history{end}(:, 2))), 'r*')
            hold off
        end
    end
    
end

