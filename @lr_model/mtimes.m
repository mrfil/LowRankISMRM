 function vo = mtimes(a, vi)
% This does the multiplication by a vector for both forward and transpose
% operations

if a.is.empty
	error empty
end

  time_basis = a.time_basis;
  ky_ind = a.ky_ind;
  N_tpt = a.N_tpt;   
  N = a.N;
  L = a.L;
       

       
   if ~a.is.transpose
       obj_pred = zeros(N*N, N_tpt);
       spat_map = reshape(vi,[N*N L]);
       for jj = 1:L
           obj_pred = obj_pred + spat_map(:,jj)*(time_basis(:,jj)');
       end
       
       % Simulation IMAGING SIGNAL == GET RANDOMIZED OR LINEAR LINE OF KSPACE
       obj_pred = reshape(obj_pred, [N N N_tpt]);
       signal_pred = ifftshift(fft(fft(fftshift(obj_pred),[],2),[],1));
       
       data_pred = zeros(N,N_tpt);
       %keyboard
       for ii = 1:N_tpt
           data_pred(:,ii) = signal_pred(ky_ind(ii),:,ii);
       end
       vo = data_pred(:);   
    
   else
       
       dat_in = reshape(vi,[N N_tpt]);    
       data_full = zeros(N, N, N_tpt);
       for ii = 1:N_tpt
           data_full(ky_ind(ii),:,ii) = dat_in(:,ii);
       end
       obj_adj = N*N*ifftshift(ifft(ifft(fftshift(data_full),[],2),[],1));
       obj_adj = reshape(obj_adj, [N*N N_tpt]);
       
       for jj = 1:L
           spat_map_up(:,jj) = obj_adj*(time_basis(:,jj));
       end
       vo = spat_map_up(:);
       
   end
