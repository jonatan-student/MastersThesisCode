classdef fdlb < fdIntegrator
	
	properties
		tau;
		cxs,cys;
		weights; 
		slipLength;
	end

	methods		
	
		function this = fdlb(tau, Ls, obstacle)
			if nargin < 1
				tau = 0.6;
			end
			this.tau = tau;
			this.cxs = [0, 0, 1, 1, 1, 0, -1, -1, -1]; 
			this.cys = [0, 1, 1, 0, -1, -1, -1, 0, 1];
			this.weights = [4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36];

			#this.slipLength = zeros(ngrdy, ngrdx);      % No slip by default
			this.slipLength = zeros(size(obstacle));        % Set slip length at walls
			this.slipLength(obstacle ==1) = Ls;
		end

		function [ux, uy, rho, f] = step(this, f, obstacle)

			[ngrdy, ngrdx, nflows] = size(f.value);
			[cidy, cidx] = find(obstacle == 1);  % solid boundary indices
			
			wall_y = [1, size(obstacle,1)];
			wall_x = 1:size(obstacle,2);

			[cidy_wall_top, cidx_wall_top] = deal(wall_y(1) * ones(1, length(wall_x)), wall_x);
			[cidy_wall_bottom, cidx_wall_bottom] = deal(wall_y(2) * ones(1, length(wall_x)), wall_x);

			% Merge with existing obstacle indices:
			cidy = [cidy; cidy_wall_top'; cidy_wall_bottom'];
			cidx = [cidx; cidx_wall_top'; cidx_wall_bottom'];
		
			% === Step 0: Apply simple boundary conditions (optional walls) ===
			f.value(:, end, [7,8,9]) = f.value(:, end-1, [7,8,9]); % Right wall
			f.value(:, 1, [3, 4, 5]) = f.value(:, 2, [3, 4, 5]);   % Left wall
			f.value(1, :, [5,6,7]) = f.value(2, :, [5,6,7]);       % Bottom wall
			f.value(end, :, [2,3,9]) = f.value(end-1, :, [2,3,9]); % Top wall
		
			% === Step 1: Streaming ===
			f.value = fdstreamlb(f.value, this.cxs, this.cys);
			f.value(1, :, :) = f.value(2, :, :);       % Bottom
			f.value(end, :, :) = f.value(end-1, :, :); % Top


			% === Step 2: Collision (compute rho, ux, uy, relax f) ===
			[f.value, rho, ux, uy] = fdcollisionlb(f.value, this.tau, this.cxs, this.cys, this.weights);
			
		
			% === Step 3: Apply bounce-back distributions ===
			[f.value, ux, uy] = fdboundaryapply(f.value, ux, uy, bndryF, cidy', cidx', this.slipLength);
		
		end
		
		
		

	end

end