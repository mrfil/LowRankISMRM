function dim = size(ob)
%function dim = size(ob)
%       "size" method for Gtomo2 class


  time_basis = ob.time_basis;
  ky_ind = ob.ky_ind;
  N_tpt = ob.N_tpt;   
  N = ob.N;
  L = ob.L;
  
  
  dim = [N*N_tpt N*N*L];