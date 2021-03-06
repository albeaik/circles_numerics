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

        % introduce the source term!!!
        %density_tau = density_tau + sourcetermvalue(triangle_i);

        %~~~~~~~~ | discretization validation and maintanance ---------------------------------------
        % current implementation of mesh and time refinemens are improper!!
        % need proper considerations for two step flow paths: step ok, step
        % redo..
        [dt_validated, DT_tau_validated, density_tau_validated, discretizationStepIsValid, deadEndFail] = solution_mesh.DiscretizationValidationAndMaintanance(dt, DT_tau, density_tau);
        
        if(deadEndFail)
            disp('Discretization scheme failed!!! solver terminated!')
            break;
        end
        
        if(discretizationStepIsValid)
            %~~~~~~~~ | step user defined coupling code -----------------------
            % Question: do we care about the precise order when this coupler code is executed? before/after remesh, etc?
            userdef_coupler.SimulateStep(dt, PDEModel, solution_mesh);
        
            %~~~~~~~~ | commit PDE step ---------------------------------------
            solution_mesh.CommitStep(dt, DT_tau_validated, density_tau_validated);
        
            %~~~~~~~~ | advance time stamp ------------------------------------
            time = time + dt;

            %~~~~~~~~ | draw --> shouldn't create side effects after this line
            if(do_realtimedraw)
                figure(fighandle)
                solution_mesh.VisualizeStep(-1);
                userdef_coupler.VisualizeStep(-1);
                drawnow
            end
        else    %Don't commit, advance time and don't draw
            dt = dt_validated;
        end
    end
    
    solution = solution_mesh;
end