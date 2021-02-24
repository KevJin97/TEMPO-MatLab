classdef bath
	properties
		dim;
		comm;
		acomm;
		WEdeg;
		NSdeg;
		eta_fun;
		eta_list;
		dt;
		Jw;
		T;
	end
	
	methods
		function obj = bath(op, Jw, T, eta)
			if nargin == 0
				obj.Jw = [];
				obj.T = 0;
                obj.dt = 0;
				obj.dim = [0, 0];
				obj.comm;
				obj.acomm;
				obj.WEdeg;
				obj.NSdeg;
			elseif Jw ~= 0
                obj.Jw = Jw;
                obj.T = T;
                obj.dt = 0;
				obj.num_eta(T, Jw)
			else
				obj.eta_fun = eta;
                obj.dt = 0;
                obj.eta_list = [];
			end
		end
		
		%function num_eta(T, Jw, subdiv)
        function num_eta(T, Jw)
			%if nargin == 2
			%	subdiv = 1000;
			%end
			function eqn = intRe(t)
				eqn = integral(@(w) w^(-2) * Jw(w) * (1 - cos(w*t)), 0, Inf);
			end
		
			function eqn = intReT(t)
				eqn = integral(@(w) w^(-2) * Jw(w) * (1 - cos(w*t)) * coth(w/(2*T)), 0, Inf);
			end
		
			function eqn = intIm(t)
				eqn = integral(@(w) w^(-2) * Jw(w) * (sin(w*t) - w*t), 0, Inf);
			end
		
			function eqn = eta(t)
				if T == 0
					eqn = intRe(t) + 1i * intIm(t);
				else
					eqn = intReT(t) + 1i * intIm(t);
                end
            end
            obj.eta_fun = eta;
		end
		
        function discretise(dt, kmax)
            if nargin ~= 2
				kmax = 0;
            end
            
            %ctime = time();
            
            if obj.dt ~= dt
                    printf('setting bath timestep')
                    obj.dt = dt;
                    obj.eta_list = [];
            end
            printf('discretising')
            %ite = list(self.eta+fun, array(range(len(self.eta_list), kmax+3))*self.dt));
            
            %for el in ite:
            %   self.eta_list.append(el);
            %print('time: ' + str(round(-ctime + time(), 2));
		end
		
        function row_degeneracy(matrix)
            %mat = array(matrix)
            mat = matrix;
            
            %v - ascontiguousarray(mat.T).view(dtype((void,
            %(mat.T).dtype.itemsize*(mat.T).shape[1])))
            %un = unique(v,return_index=True, return_inverse=True)
            %return [len(un[0]), array(un[1], un[2])
		end
		
        function I_dk(dk, unique)
            
		end
	end
end
