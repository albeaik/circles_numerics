% This file contains pre-conditions used to test and evaluate our solvers


Q0_1 = @(x,v) (x >= 10) .* (x <= 30) .* (v >= 5) .* (v <= 10);
Q0_2 = @(x,v) (x>=2).*(x<=5).*(v>=10).*(v<=15) + (x>=5).*(x<=15).*(v>=5).*(v<=10);
Q0_3 = @(x,v) (x >= 8) .* (x <= 12) .* (v >= 7) .* (v <= 12);
Q0_4 = @(x,v) (x>=2).*(x<=5).*(v>=12).*(v<=13) + (x>=5).*(x<=15).*(v>=5).*(v<=7);
Q0_5 = @(x,v) (x>=0).*(x<=10).*(v>=12).*(v<=13) + (x>=10).*(x<=20).*(v>=5).*(v<=7);
Q0_6 = @(x,v) (x >= 8) .* (x <= 12);



[Q0_7.dbp1, Q0_7.dbe1] = get_rectangler_discontinuity_boundary( [0,12], [10,13], [0.5,1] ); 
[Q0_7.dbp2, Q0_7.dbe2] = get_rectangler_discontinuity_boundary( [10,5], [20,7], [0.5,1] );
Q0_7.oracle = @(x,v) (x>=0).*(x<=10).*(v>=12).*(v<=13) + (x>=10).*(x<=20).*(v>=5).*(v<=7);
Q0_7.discontinuity_boundary_points = [Q0_7.dbp1; Q0_7.dbp2];
Q0_7.discontinuity_boundary_edges = [Q0_7.dbe1; Q0_7.dbe2+size(Q0_7.dbp1, 1)];

[Q0_8.dbp1, Q0_8.dbe1] = get_rectangler_discontinuity_boundary( [0,12], [10,13], [0.5,1] ); 
Q0_8.oracle = @(x,v) (x>=0).*(x<=10).*(v>=12).*(v<=13);
Q0_8.discontinuity_boundary_points = [Q0_8.dbp1];
Q0_8.discontinuity_boundary_edges = [Q0_8.dbe1];