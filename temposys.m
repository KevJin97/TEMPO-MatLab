classdef temposys
   properties
       dim;
       ham;
       b;
       state;
       istate;
       freeprop;
       dmax;
       prec;
       point;
       dt;
       statedat;
       name;
   end
   
   methods
       function obj = temposys(hilbert_dim)
           if nargin == 0
               hilbert_dim = 2;
           end
           
           obj.dim = hilbert_dim;
           obj.ham = zeroes(org.dim^2, org.dim^2);
           obj.b = [];
           obj.state = zeroes(obj.dim^2);
           obj.istate = zeroes(obj.dim^2);
           obj.freeprop = zeroes(obj.dim^2, obj.dim^2);
           obj.dkmax = 0;
           obj.prec = 0;
           obj.dt = 0;
           obj.point = 0;
           obj.statedat = [ [], [] ];
           obj.name = 'temp';
       end
       
       function superop = comm(op)
           superop;
           %return kron(op, eye(self.dim)) - kron(eye(self.dim), op.T)
           %%%%%transposed
       end
       
       function superop = acomm(op)
           
       end
       
       function superop = diss(op)
           
       end
       
       function set_filename(name_string)
           
       end
       
       function checkdim(op_array)
           
       end
       
       function set_state(state_array)
           
       end
       
       function set_hamiltonian(ham)
           
       end
       
       function add_dissipation(gam, lind)
           
       end
       
       function get_state()
           
       end
       
       function convergence_params(dt, dk, prec)
           
       end
       
       function add_bath(op, Jw, T, eta)
           
       end
       
       function data = getopdat(op)
           
       end
       
       
   end
end