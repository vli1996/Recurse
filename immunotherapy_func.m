function [posA, posB, posC, livedeadB, livedeadC, livedeadsum, t ] = immunotherapy_func( n0_A, n0_B, n0_C, speedA, speedB, dt, tot_steps, d_crit_sq, kd, growth, evolve, dims )
% A is the immune cells, that have speedA
% B is the non-evolved cancer cells, that have speedB
% C is the evolved cancer cells, that have speedB
% kd is the rate at which immune cells destroy tumor cells
% d_crit_sq is the squared distance between immune cells and tumor cells to
% initiate reaction
% growth is the rate at which tumor cells multiply
% evolve is the rate at which tumor cells evolve

t = 0:dt:dt*tot_steps;

posA = zeros(length(t), n0_A, 2); % keeps track of the position of all, rows r time steps, columns r molecules
posB = zeros(length(t), n0_B, 2); 
posC = zeros(length(t), n0_C, 2); 
livedeadB = ones(length(t), n0_B); % keeps track of which molecules (columns) r alive at what times (rows)
livedeadC = ones(length(t), n0_C);
livedeadsum = zeros(1, length(t));
livedeadsum(1) = n0_B + n0_C;

% random placement of cells

for i = 1:n0_A
    posA(1, i,:) = 100*[rand(1) rand(1)] - 50;
end

for i = 1:n0_B
    posB(1, i,:) = 100*[rand(1) rand(1)] - 50;
end

for i = 1:n0_C
    posC(1, i,:) = 100*[rand(1) rand(1)] - 50;
end

targetarray = randperm(n0_B);
targetindex = targetarray(1:n0_A); % select random cancer cells to target, one for each immune cell

% movement
for i = 1:length(t) - 1
    
    % loop sizes
    loopsizeA = size(posA);
    loopsizeB = size(posB);
    loopsizeC = size(posC);
    
    % movement of immune cells
    
    for j = 1:loopsizeA(2)
        for k = 1:dims
            targetpos = posB(i, targetindex(j), k);
            
            if posA(i, j, k) < targetpos
                posA(i+1, j, k) = posA(i, j, k) + (speedA*dt)^0.5;
            else
                posA(i+1, j, k) = posA(i, j, k) - (speedA*dt)^0.5;
            end
            % posA(i+1, j, k) = posA(i,j,k) + (speedA*dt)^0.5*randn(1); % da fuq + (speedA*dt) * posB(i, ; % needs movement component 
        end
    end
    
    
    
    % movement of tumor cells
    for j2 = 1:loopsizeB(2)
        
        % new tumor cell growth, this doesn't work yet
        if rand(1) < growth*dt & livedeadB(i, j2) ~= 0
            newB = zeros(length(t), 1, 2);
            newB(i + 1, 1, :) = posB(i, j2, :); % sets position of the new tumor cell to = the old one, doesn't work
            posB = [posB newB]; % adds a new cell to keep track of in posB
            newBlivedead = zeros(length(t), 1);
            newBlivedead(i+1:end, 1) = 1;
            livedeadB = [livedeadB newBlivedead]; % adds a new column with 0's leading up to the growth
        end
        
        % evolution
        if rand(1) < evolve*dt & livedeadB(i, j2) ~= 0
            newC = zeros(length(t), 1, 2);
            newC(i + 1, 1, :) = posB(i, j2, :); 
            posC = [posC newC]; 
            livedeadB(i+1:end,j2) = 0;
            newClivedead = zeros(length(t), 1);
            newClivedead(i + 1:end, 1) = 1;
            livedeadC = [livedeadC newClivedead];
        end
        
        for k = 1:dims
            posB(i+1, j2, k) = posB(i,j2,k) + (speedB*dt)^0.5*randn(1);
        end
        
        % immune cell kill
        for j3 = 1:n0_A
            if ((posB(i+1,j2,1) - posA(i+1,j3,1))^2 + (posB(i+1,j2,2) - posA(i+1,j3,2))^2) <= d_crit_sq
                if rand(1) < kd*dt  
                    livedeadB(i+1:end,j2) = 0;
                end
            end
            
            if livedeadB(i+1, targetindex(j3)) == 0
                    randarray = randperm(loopsizeB(2));
                    targetindex(j3) = randarray(1);
            end
        end
    end
    
    % evolved tumor cell movement
    for j4 = 1:loopsizeC(2)
        % evolved tumor cell growth
        if rand(1) < growth*dt & livedeadC(i, j4) ~= 0
            newC = zeros(length(t), 1, 2);
            newC(i + 1, 1, :) = posC(i, j4, :);
            posC = [posC newC];
            newClivedead = zeros(length(t), 1);
            newClivedead(i+1:end, 1) = 1;
            livedeadC = [livedeadC newClivedead];
        end
        
        for k = 1:dims
            posC(i+1, j4, k) = posC(i,j4,k) + (speedB*dt)^0.5*randn(1);
        end  
    end
    
    livedeadsum(i+1) = sum(livedeadB(i+1,:)) + sum(livedeadC(i+1, :));
    
end

end