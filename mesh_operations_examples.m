close all
clear all
clc


load('../simulations/july_2/sample')
[ DT_new ] = ExpandMeshedDomain( q_solution_mesh.DT, mesh_domain_limits, mesh_resolution ); 
figure
subplot(2, 1, 1)
trimesh(q_solution_mesh.DT.ConnectivityList, q_solution_mesh.DT.Points(:, 1), q_solution_mesh.DT.Points(:, 2))
subplot(2, 1, 2)
trimesh(DT_new.ConnectivityList, DT_new.Points(:, 1), DT_new.Points(:, 2))


pre_specified_initial_conditions;

mesh_domain_limits = [0, 20; 0, 20];
mesh_resolution = [40; 40];
mesh_1 = DiscreteTimeEvolvingMesh(Q0_3, mesh_domain_limits, mesh_resolution);

mesh_domain_limits = [0, 20; 0, 20];
mesh_resolution = [45; 45];
mesh_2 = DiscreteTimeEvolvingMesh(Q0_2, mesh_domain_limits, mesh_resolution);


discont_1 = get_discontinuity_edge_list(mesh_1.DT, mesh_1.density, 0.5);
discont_2 = get_discontinuity_edge_list(mesh_2.DT, mesh_2.density, 0.5);

new_mesh_points = [mesh_1.DT.Points; mesh_2.DT.Points];
discont_1;
discont_2 = discont_2 + size(mesh_1.DT.Points, 1);
new_discont = [discont_1; discont_2];

new_DT = delaunayTriangulation(new_mesh_points, new_discont);
density_contribution_mesh_1 = mesh_1.density(mesh_1.DT.pointLocation(new_DT.incenter));
density_contribution_mesh_2 = mesh_2.density(mesh_2.DT.pointLocation(new_DT.incenter));
new_density = density_contribution_mesh_1 + density_contribution_mesh_2;

figure
subplot(3, 1, 1)
DiscreteTimeEvolvingMesh.visualize_triangulation(mesh_1.DT, mesh_1.density);
caxis([0, 2]);
subplot(3, 1, 2)
DiscreteTimeEvolvingMesh.visualize_triangulation(mesh_2.DT, mesh_2.density);
caxis([0, 2]);
subplot(3, 1, 3)
DiscreteTimeEvolvingMesh.visualize_triangulation(new_DT, new_density);
caxis([0, 2]);
