close all
clear all
clc



nT = 100;                %number of time steps
nAC = 2;                %number of autonomous cars
u0 = zeros([nT*nAC, 1]);

fun = @(u) optimal_control_simple_obj_fcn(u, nAC, nT);

min_accel = 0;
max_accel = 6;

A = [-1*eye(nT*nAC); eye(nT*nAC)];
b = [-min_accel*ones([nT*nAC, 1]); max_accel*ones([nT*nAC, 1])];

u_opt = fmincon(fun, u0, A, b);