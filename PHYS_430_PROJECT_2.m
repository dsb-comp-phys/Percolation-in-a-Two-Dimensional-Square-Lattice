clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EDIT THESE VARIABLES TO EXPERIMENT WITH THE PROGRAM!

L = 25;             % SET SIZE OF LATTICE HERE (L>300 not reccomended)
ANIMATE = true;     % BOOLEAN TO SHOW ANIMATION (L>40 not reccomended)

%{
This program makes use of additional MATLAB packages. If you do not have
these packages installed an run into a related error with the program, set
the following variable to false. Otherwise, leave this as true.
%}
NON_LINEAR_FIT = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

N = L^2;
N_rem = N;

rem_array = 1:N; % remaining pos to be filled

lattice = zeros(L);

% At first, F_c counts the number of cells in the percolating cluster.
% Only at the end will it be divided by N, whence being a fraction.
F_c_array = zeros(1,N);

% This is the probability at which a percolating cluster occurs.
% Note that when we speak of probability, we speak of cells/N.
% Cells, meaning filled cells.
p_c = 0;
p_c_cell = 0; % cell at which p_c is found

perc_found = false;

% Boolean to determine whether the plot for p=0.8 has already been created. 
checked80 = false;

% Loop of cells
for cells = 1:N
    
    if cells ~= 1
        F_c_array(cells) = F_c_array(cells-1);
    end
	
	r = randi([1,N_rem]); % position in rem_array
	pos = rem_array(r);   % position in lattice
	rem_array(r) = [];    % remove element from rem_array
	N_rem = N_rem - 1;    % length of rem_array has decreased by 1
	
	% Note that we are working with a square lattice.
	% Hence, we can get the row and column from pos.
	[i,j] = pos_ij(pos,L);
	
	% Finally, we fill the cell found.
	lattice(i,j)=1;
	
	% 0: UNOCCUPIED CELL.
	% 1: OCCUPIDED CELL, BUT NOT PART OF PERCOLATING CLUSTER.
	% 2: OCCUPIED CELL THAT IS PART OF THE PERCOLATING CLUSTER.
	
	% If a percolating cluster has yet to be found...
	% From the new position, we check for a percolating cluster.
	% We use our BFS function for this.
	% If so, perc_pos = pos;
	% Save the current p as the percolating probability p_c.
	
	% If a percolating cluster has been found in this iteration...
	% We start from per_pos and count the number of cells in its cluster.
	% To do this, we perform BFS.
	% All cells in the percolating clusters will be given the value 2.
	% 2 is assigned when the cell is marked as visited.
	% Use the to find the fraction of cells F_c in percolating cluster.
	
	% If a percolating cluster was found in a previous iteration...
	% If the current pos has a neighbor with value 2, run BFS.
	% Add newly visted cells to the count of cells in percolating cluster.
	% Update F_c.
	
    run_bfs = false;
    
    if perc_found
        
        [up,left,right,down] = neighbors(pos,L);
        
        [up_i,up_j] = pos_ij(up,L);
        [left_i,left_j] = pos_ij(left,L);
        [right_i,right_j] = pos_ij(right,L);
        [down_i,down_j] = pos_ij(down,L);
        
        if up~=0
            if lattice(up_i,up_j)==2
                run_bfs = true;
            end
        end
        
        if left~=0
            if lattice(left_i,left_j)==2
                run_bfs = true;
            end
        end
        
        if right~=0
            if lattice(right_i,right_j)==2
                run_bfs = true;
            end
        end
        
        if down~=0
            if lattice(down_i,down_j)==2
                run_bfs = true;
            end
        end
        
    else
        run_bfs = true;     
    end
    
    if run_bfs
        [lattice_update, F_c_update] = bfs(pos, lattice, F_c_array(cells));
        
        if F_c_update ~= 0 && ~perc_found
            perc_found = true;
            p_c_cell = cells;
            p_c = cells/N;
            F_c_array(cells) = F_c_update;
            F_cell = cells;
            lattice = lattice_update;
            
            figure(2)
            imagesc(lattice)
            h = colorbar;
            set(h, 'ylim', [0 2])
            title(['p_c=', num2str(round((N-N_rem)/N,3)),' , F_c=', num2str(round(F_c_array(cells)/cells,3))])
            xlabel('x')
            ylabel('y')
            pbaspect([1 1 1])
        end
    
        if perc_found
            F_c_array(cells) = F_c_update;
            lattice = lattice_update;
        end
        
    end
	
	% SURFACE PLOT OF LATTICE
    if ANIMATE
        figure(1)
        imagesc(lattice)
        h = colorbar;
        set(h, 'ylim', [0 2])
        title(['p=', num2str(round((N-N_rem)/N,3)),' , F_c=', num2str(round(F_c_array(cells)/cells,3))])
        xlabel('x')
        ylabel('y')
        pbaspect([1 1 1])
    end
    
    if round(cells/N,2) == 0.80
        if ~checked80
            
            checked80 = true;
            
            figure(3)
            imagesc(lattice)
            h = colorbar;
            set(h, 'ylim', [0 2])
            title(['p=', num2str(round((N-N_rem)/N,3)),' , F_c=', num2str(round(F_c_array(cells)/cells,3))])
            xlabel('x')
            ylabel('y')
            pbaspect([1 1 1])
        end
    end
	
end     % cells loop

% THIS IS A VERY IMPORTANT LINE OF CODE!
% We are finding the fraction of percolating cells to occupied cells
% Hence, (cells labeled 2) / ( (cells labeled 2) + (cells labeled 1) )
F_c_array = F_c_array ./ (1:N);

p=(1:N)/N;

disp("p_c="+p_c)

fun = @(F_0,beta,x) F_0*(x-p_c).^beta;

if NON_LINEAR_FIT
    fit_plot=fit(p(1,p_c_cell:N)',F_c_array(1,p_c_cell:N)',fun,'Lower',[.1,0])
end


figure(4)
plot((1:N)/N,F_c_array, '.')
hold on
%plot((1:N)/N,F_0*((1:N)/N-p_c).^beta)
if NON_LINEAR_FIT
    plot(fit_plot)
end
hold off
title(['2D percolation ', num2str(L), '\times', num2str(L), ' lattice'])
xlabel('p')
ylabel('F_c')
xlim([p_c 1])
ylim([0 1])



% i and j are found for the given pos in the sqare lattice.
function [i,j] = pos_ij(pos,L)

    i = ceil(pos/L);
    j = rem(pos,L);

    if j == 0
        j = L;
    end

end


function pos = ij_pos(i,j,L)

    pos = (i-1)*L+j;
    
end


function [up,left,right,down] = neighbors(pos,L)

    [i,j] = pos_ij(pos,L);
    
    if i~=1
        up = ij_pos(i-1,j,L);
    else
        up = 0;
    end
    
    if j~=1
        left = ij_pos(i,j-1,L);
    else
        left = 0;
    end
    
    if j~=L
        right = ij_pos(i,j+1,L);
    else
        right = 0;
    end
    
    if i~=L
        down = ij_pos(i+1,j,L);
    else
        down = 0;
    end
    
end


function [lattice_update, F_c_update] = bfs(pos, lattice, F_c)
    
    % Booleans are used to find if the cluster percolates.
    top_row = false;
    left_column = false;
    right_column = false;
    bottom_row = false;
    
    percolates = false;

    lattice_update = lattice;
    F_c_update = F_c;
    
    L = length(lattice_update);

    % THIS IS PYTHONIC, BUT I'M DOING THIS ANYWAY.    
    q = [pos];
    
    [pos_i,pos_j] = pos_ij(pos,L);
    
    lattice_update(pos_i,pos_j) = 2;
    F_c_update = F_c_update + 1;

    % This algorithm prioritizes visiting neighbors in the following way:
    % up, left, right, down (by value of respective pos).

    while ~isempty(q)
        
        % NOTE: perc_found is a boolean for whether a percolating cluster 
        % has already been found.
        
        % Define new pos
        pos = q(1);
        q(1) = [];  % dequeue new pos
        
        % Define pos neighbors
        % Label 0 if non-existent
        [up,left,right,down] = neighbors(pos,L);
        
        [up_i,up_j] = pos_ij(up,L);
        [left_i,left_j] = pos_ij(left,L);
        [right_i,right_j] = pos_ij(right,L);
        [down_i,down_j] = pos_ij(down,L);
        
        if up == 0
            top_row = true; % cluster inhabits top row
        elseif lattice_update(up_i,up_j) == 1
                lattice_update(up_i,up_j) = 2;
                F_c_update = F_c_update + 1;
                q(length(q)+1)=up;
        end
    
        if left == 0
            left_column = true; % cluster inhabits leftmost column
        elseif lattice_update(left_i,left_j) == 1
                lattice_update(left_i,left_j) = 2;
                F_c_update = F_c_update + 1;
                q(length(q)+1)=left;
        end
        
        if right == 0
            right_column = true; % cluster inhabits rightmost column
        elseif lattice_update(right_i,right_j) == 1
                lattice_update(right_i,right_j) = 2;
                F_c_update = F_c_update + 1;
                q(length(q)+1)=right;
        end
        
        if down == 0
            bottom_row = true; % cluster inhabits bottom row
        elseif lattice_update(down_i,down_j) == 1
                lattice_update(down_i,down_j) = 2;
                F_c_update = F_c_update + 1;
                q(length(q)+1)=down;
        end
                
        % Enqueue neighbors and mark as visited
            
    end % while
    
    if (top_row && bottom_row) || (left_column && right_column)
        percolates = true;
    end
    
    if ~percolates && F_c == 0
        % A percolating cluster has yet to be found
        lattice_update = lattice; % no update to lattice afterall
        F_c_update = F_c;         % nor to F_c
    end

end % bfs
