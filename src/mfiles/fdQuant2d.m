
classdef fdQuant2d < fdQuantity 

	properties (Access=public)
		# Grid dimensions
		ngrdx, ngrdy; 
		
		# Grid spacings
		dx, dy;
		
		# Derivative
		ddx; ddy;
		d2dx2, d2dy2;
		bcvals;
	end


	methods  

		function this = fdQuant2d(grids, spacings, boundaries)
			
			if ( nargin == 1 )
				this;
			elseif ( nargin >= 2 )

				if ( length(grids) != 2 || length(spacings) != 2 )
					error("Specification of the grid and spacings not valid");
				end

				this.ngrdx = grids(1); this.ngrdy = grids(2);
				this.dx = spacings(1); this.dy = spacings(2);
				this.value = this.nvalue = this.rhs = zeros(this.ngrdy, this.ngrdx);
				
			end
			
			if ( nargin == 3 )
			
				if ( isstring(boundaries) != 0 || length(boundaries) != 4 )
					error("Boundary specification not valid");
				end

				this.bcs = boundaries;
			end	
		end

		function lbset(this, grids)
			this.ngrdx = grids(1); this.ngrdy = grids(2);
			this.value = ones(this.ngrdy, this.ngrdx, 9)+ 0.001*randn(this.ngrdy, this.ngrdx);
		end


		function retval = calcddx(this, method)
			
			if nargin==1 
				this.ddx = fdcds2d(this.value, this.dx, 'x');
			else 
				switch(method)
					case "central"
						this.ddx = fdcds2d(this.value, this.dx, 'x');
					case "forward"
						this.ddx = fdfrwd2d(this.value, this.dx, 'x');
					case "backward"
						this.ddx = fdbckwd2d(this.value, this.dx, 'x');
					otherwise
						error("Invalid method");
				end 
			end

			if this.bcs(2)=='p' 
				this.ddx(:, 1) = (this.value(:,2) - this.value(:,end))./(2*this.dx);
    			this.ddx(:, end) = (this.value(:,1) - this.value(:,end-1))./(2*this.dx);
			end

			retval = this.ddx;
		end

		function retval = calcddy(this, method)
			
			if nargin==1
				this.ddy = fdcds2d(this.value, this.dy, 'y');
			else			
				switch(method)
					case "central"
						this.ddy = fdcds2d(this.value, this.dy, 'y');
					case "forward"
						this.ddy = fdfrwd2d(this.value, this.dy, 'y');
					case "backward"
						this.ddy = fdbckwd2d(this.value, this.dy, 'y');
					otherwise
						error("Invalid method");	
				end 
			end
			
			if this.bcs(1)=='p'
				this.ddy(1,:) = (this.value(2,:) - this.value(end,:))./(2*this.dy);
    			this.ddy(end,:) = (this.value(1,:) - this.value(end-1,:))./(2*this.dy);
			end

			retval = this.ddy;
		end

		function retval = calcd2dx2(this)
			
			this.d2dx2 = fdscnd2d(this.value, this.dx, 'x');

			if this.bcs(2)=='p'
				idx2 = 1.0/(this.dx^2);
				this.d2dx2(:, 1) = (this.value(:,2) - 2*this.value(:,1) + this.value(:,end))*idx2;
			    this.d2dx2(:, end) = (this.value(:,1) - 2*this.value(:,end) + this.value(:,end-1))*idx2;
			end 

			retval = this.d2dx2;
		end

		function retval = calcd2dy2(this)
			
			this.d2dy2 = fdscnd2d(this.value, this.dy, 'y');

			if this.bcs(1)=='p'
				idy2 = 1.0/(this.dy^2);
				this.d2dy2(1,:) = (this.value(2,:) - 2*this.value(1,:) + this.value(end,:))*idy2;
			    this.d2dy2(end,:) = (this.value(1,:) - 2*this.value(end,:) + this.value(end-1,:))*idy2;
			end 
			
			retval = this.d2dy2;
		end

		function retval = laplace(this)
			
			retval = this.calcd2dx2() + this.calcd2dy2();			

		end

		function L = laplaceMatrix(this)
			ny = this.ngrdy;
			nx = this.ngrdx;
			N = ny * nx;
		
			dx2 = this.dx^2;
			dy2 = this.dy^2;
		
			Ix = speye(nx);
			Iy = speye(ny);
		
			e = ones(nx,1);
			Tx = spdiags([e -2*e e], -1:1, nx, nx) / dx2;
			Ty = spdiags([e -2*e e], -1:1, ny, ny) / dy2;
		
			L = kron(Ix, Ty) + kron(Tx, Iy);  % 2D Laplacian
		end
		

		function setvalue(this, value)
			this.value = value*ones(this.ngrdy, this.ngrdx);
			this.value(1,:) = this.value(this.ngrdy, :) = 0; 
			this.value(:,1) = this.value(:,this.ngrdx) = 0;
		end

		function setbcs(this, specifier)
			if ( !ischar(specifier) || length(specifier) != 4 )
				error("Boundary specifications not correct");
			end
			this.bcs = specifier;
		end

		function setbc(this, side, type, value)
			% side  : 'top'|'right'|'bottom'|'left'
			% type  : 'dirichlet'|'neumann'|'periodic'
			% value : scalar (used only when type == 'dirichlet')
		
			names = {'top','right','bottom','left'};
			idx   = find(strcmpi(side,names));
			if isempty(idx),  error('fdQuant2d:setbc','Bad side string');  end
		
			switch lower(type)
				case 'dirichlet', this.bcs(idx) = 'd';
				case 'neumann',   this.bcs(idx) = 'n';
				case 'periodic',  this.bcs(idx) = 'p';
				otherwise,        error('fdQuant2d:setbc','Bad BC type');
			end
		
			if isempty(this.bcvals),  this.bcvals = zeros(1,4);  end
			this.bcvals(idx) = value;
		end
		

		function applybcs(this)

			% make sure bcvals exists
			if isempty(this.bcvals),  this.bcvals = zeros(1,4);  end
		
			for n = 1:4
				switch this.bcs(n)
		
					% -------- Dirichlet -----------------------------------------
					case 'd'
						val = this.bcvals(n);
						if     n == 1      % top
							this.value(1,:)            = val;
						elseif n == 2      % right
							this.value(:, this.ngrdx)  = val;
						elseif n == 3      % bottom
							this.value(this.ngrdy,:)   = val;
						else               % left
							this.value(:,1)            = val;
						end
		
					% -------- Neumann (zero-gradient) ---------------------------
					case 'n'
						if     n == 1
							this.value(1,:)           = this.value(2,:);
						elseif n == 2
							this.value(:, this.ngrdx) = this.value(:, this.ngrdx-1);
						elseif n == 3
							this.value(this.ngrdy,:)  = this.value(this.ngrdy-1,:);
						else
							this.value(:,1)           = this.value(:,2);
						end
		
					% -------- Periodic: nothing to enforce here -----------------
					case 'p'
						% derivatives already wrap around; leave values as-is
				end
			end
		end
		

		function update(this)
			this.pvalue = this.value;
			this.value = this.nvalue;		

			this.applybcs();	
		end


	end


end
