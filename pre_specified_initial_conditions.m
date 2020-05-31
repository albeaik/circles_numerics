% This file contains pre-conditions used to test and evaluate our solvers


Q0_1 = @(x,v) (x >= 10) .* (x <= 30) .* (v >= 5) .* (v <= 10);
Q0_2 = @(x,v) (x>=2).*(x<=5).*(v>=10).*(v<=15) + (x>=5).*(x<=15).*(v>=5).*(v<=10);