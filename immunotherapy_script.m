
% Immune cells are displayed in red and track tumor cells, destroying them
% in a certain proximity

% Tumor cells move through a random walk, can multiply, and can also evolve
% into resistant tumor cells which can't be tracked or destroyed by immune
% cells. These evolved immune cells also move and multiply

clear all;
close all;

speedA = 2;
speedB = 1;
dt = 1;
n0_A = 5;
n0_B = 25;
n0_C = 0;

tot_steps = 100;
dims = 2;

d_crit_sq = 3^2;
kd = .25;
growth = 0.005;
evolve = 0.005;

[posA, posB, posC, livedeadB, livedeadC, livedeadsum, t ] = immunotherapy_func( n0_A, n0_B, n0_C, speedA, speedB, dt, tot_steps, d_crit_sq, kd, growth, evolve, dims );

figure;
plot(t, livedeadsum);
title('Number of Cancer Cells Over Time');
axis([0 dt*tot_steps 0 n0_B]);


figure;
for i = 1:length(t)
    hold off;
    for j = 1:n0_A
        plot(posA(i,j,1), posA(i,j,2), 'o', 'color', 'r'); % plots immune cell position
        hold on;
    end
    
    title('Immunotherapy Simulation');
    xlabel('Red = Immune Cells, Blue = Unevolved Tumor Cells, Green = Evolved Tumor Cells');
    
    loopsize = size(posB);
    for j2 = 1:loopsize(2) - 1
        if livedeadB(i,j2) ~= 0 % only displays cancer cells when the time step value in the alive column is 1
            plot(posB(i,j2,1), posB(i,j2,2), 'o', 'color', 'b');
        end
    end
    
    loopsizeC = size(posC);
    for j3 = 1:loopsizeC(2)
        if livedeadC(i,j3) ~= 0
            plot(posC(i,j3,1), posC(i,j3,2), 'o', 'color', 'g'); % plots evolved tumor cell position
        end
    end
    
    xlim([-120 120]);
    ylim([-120 120]);
    pause(.1)
end
            
        
        




