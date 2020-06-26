function [solution] = characteristic_solver(solution_obj, dt, T, G_fcn, H_fcn, user_defined_coupler_obj, do_realtimedraw)
    
    %--- if not initialized, then ??
    solution_mesh = solution_obj;
    userdef_coupler = user_defined_coupler_obj;

    %--- solve
    time = 0;       %time: seconds
    while(time <= T)
    % write down iteration contract: assume DT_t, density_t set properly, never changed until end of iteration, etc

        %~~~~~~~~ | prepare for this iteration
        time
        DT_t = solution_mesh.DT;
        density_t = solution_mesh.density;
        
        %~~~~~~~~ | step PDE: solve PDE evolution in time -----------------
        % simulate the characteristic equation
        DT_tau_Points = DT_t.Points + dt*[G_fcn(solution_mesh, userdef_coupler, DT_t.Points), ...
                                            H_fcn(solution_mesh, userdef_coupler, DT_t.Points)];

        %create DT_tau object ==> note matlab native DT object recomputes triangulation if points updated directly                       
        DT_tau.Points = DT_tau_Points;
        DT_tau.ConnectivityList = DT_t.ConnectivityList;

        % simulate density scaling
        [DT_t_trig_areas] = GetDelaunayAreas(DT_t);  %sorted by triangle id --> can pre-compute if DT is fixed
        [DT_tau_trig_areas] = GetDelaunayAreas(DT_tau);  %sorted by triangle id
        density_tau = (DT_t_trig_areas ./ DT_tau_trig_areas) .* density_t;

        %~~~~~~~~ | step user defined coupling code -----------------------
        userdef_coupler.SimulateStep(dt, solution_mesh);

        %~~~~~~~~ | commit PDE step ---------------------------------------
        solution_mesh.EvolveMesh(dt, DT_tau, density_tau);
        
        %~~~~~~~~ | advance time stamp ------------------------------------
        time = time + dt;
        
        %~~~~~~~~ | draw --> shouldn't create side effects after this line
        if(do_realtimedraw)
            solution_mesh.VisualizeMesh();
            userdef_coupler.VisualizeStep();
            drawnow
        end
    end
    
    solution = solution_mesh;
end