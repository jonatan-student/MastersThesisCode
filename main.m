clear; clc;
addpath("../../src/mfiles");
addpath("../../src/cppfiles");
pkg load parallel;

% decide how many workers
keepfree = 2;
nc = nproc();
numworkers = max(nc - keepfree, 1)
fprintf("\nTotal Processors: %d\nIn use: %d\n", nc, numworkers)

%set up parameters for tests
height_mult = [1,3,5,7,9];
Ls_array = [0, 4, 8, 12, 16, 20];
slopes_array = [0, 2, 4, 6, 8, 10]; % multiplyier by 0.01 on both top and bottom, ie. 2x0.01 = 0.02 X 2(walls) = 0.04=m
output_dir = "AfterHandinTests/toHalfPeDaTrue";
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

runMembraneRefactored(0, 'Slopes', output_dir);

%dispatch in parallel stopped working for me.
%pararrayfun(numworkers, @(var) runMembraneRefactored(var, "SlipLengths", output_dir),  Ls_array);

fprintf("\n-----------------\nSimulation complete!\n-----------------\n");
