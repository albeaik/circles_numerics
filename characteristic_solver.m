function [solution] = characteristic_solver(solution_obj, dt, T, PDEModel, user_defined_coupler_obj, do_realtimedraw)
    
    %--- if not initialized, then ??
    solution_mesh = solution_obj;
    userdef_coupler = user_defined_coupler_obj;
    
    if(do_realtimedraw)
        fighandle = figure;
    end

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
        DT_tau_Points = DT_t.Points + dt*[PDEModel.G(solution_mesh, userdef_coupler, DT_t.Points), ...
                                            PDEModel.H(solution_mesh, userdef_coupler, DT_t.Points)];

        %create DT_tau object ==> note matlab native DT object recomputes triangulation if points updated directly                       
        DT_tau.Points = DT_tau_Points;
        DT_tau.ConnectivityList = DT_t.ConnectivityList;

        % simulate density scaling
        density_tau = (solution_mesh.GetMeshAreas(DT_t) ./ solution_mesh.GetMeshAreas(DT_tau)) .* density_t;

        %~~~~~~~~ | step user defined coupling code -----------------------
        userdef_coupler.SimulateStep(dt, PDEModel, solution_mesh);

        %~~~~~~~~ | commit PDE step ---------------------------------------
        [newdt, newMesh, newMeshValues, forceUndoStep, terminate] = solution_mesh.EvolveMesh(dt, DT_tau, density_tau);
        
        if(terminate)
            disp('Discretization scheme failed!!! solver terminated!')
            break;
        end
        
        if(forceUndoStep)
            dt = newdt;
            continue;   %Don't advance time and don't drawing
        end
        
        %~~~~~~~~ | advance time stamp ------------------------------------
        time = time + dt;
        
        %~~~~~~~~ | draw --> shouldn't create side effects after this line
        if(do_realtimedraw)
            figure(fighandle)
            solution_mesh.VisualizeStep(-1);
            userdef_coupler.VisualizeStep(-1);
            drawnow
        end
    end
    
    solution = solution_mesh;
end