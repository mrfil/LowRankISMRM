%% Low Rank Reconstruction Tutorial for ISMRM 2018
% Presented by Brad Sutton, June 17, 2018
%
% This code creates a dynamic object with time varying frequencies 
% in the data and then simulates navigator and imaging data streams. 
% It then performs both a sliding window reconstruction to show the 
% severe undersampling and a low rank reconstruction that recovers 
% the spatiotemporal dynamics. This code is for illustrative purposes
% only and does not represent the best way to do low rank reconstructions
% which should include spatialspectral sparsity or other penalties to 
% accomodate the more complex real-imaging data scenarios.
%
% Code is included from Jeff Fessler's IRT package which can be found at:
% https://web.eecs.umich.edu/~fessler/code/

addpath('Fessler_IRT/')
flag_make_movies = 1;  % This determines if the movies will be generated


N = 64;  % Image matrix size

dT = 0.1;  % 10 frames per second, just an arbitrary number for this simulation
N_frames = 20;  % Number of full frames of data that will be collected
				%    20 full frames results in 20*64 lines of imaging data acquired
N_tpt = N*N_frames   % Number of time points total in the time series


%% Object simulation
% Time-varying frequencies of motion for the inner and outer objects
f_in = 2*ones(N_tpt,1)+0.3*sin(2*pi*[1:N_tpt]'/100);   
f_out = 0.3*ones(N_tpt,1)+0.15*cos(2*pi*[1:N_tpt]'/100);

% Location of centers of ovals for inner and outer objects
cx_in = -4;
cy_in = 0;

% Setting up spatial coordinates for object simulation
[y,x] = ndgrid([1:N]-N/2);


% Integrating frequencies to get phases of motion for the inner and outer objects
phi_in = cumtrapz(f_in)*dT;
phi_out = cumtrapz(f_out)*dT;

% Simulate object at each of the sampling locations in time
for ii = 1:N_tpt
    
    r_x_out = 27*ones(N_tpt,1)+4*sin(phi_out);
    r_y_out = 20*ones(N_tpt,1)+3*sin(phi_out);
    r_x_in = 14*ones(N_tpt,1)+2*sin(phi_in);
    r_y_in = 12*ones(N_tpt,1)+3*sin(phi_in);
       
    % is inside:
    mask_outside = ( ((x.^2./r_x_out(ii)^2)+(y.^2./r_y_out(ii)^2)) <1);
    mask_inside = ( (((x-cx_in).^2./r_x_in(ii)^2)+((y-cy_in).^2./r_y_in(ii)^2)) <1);
   
    obj(:,:,ii) = xor(mask_outside,mask_inside);
    
end

obj = double(obj);

%%% ADD NOISE
%obj = obj + 0.1*randn(N,N,N_tpt);

figure
if flag_make_movies    % WRITE OUT MOVIE
H = VideoWriter('Orig_Obj','MPEG-4');
H.FrameRate=10;
open(H)
for ii = 1:N_tpt
    imagesc(obj(:,:,ii),[0 1.25])
    colormap gray
    %pause(0.1)
    currframe=getframe;
    writeVideo(H,currframe);
end
close(H)
end


%% PERFORM SVD ON ENTIRE GROUND TRUTH IMAGE
X_obj = reshape(obj,[N*N N_tpt]);
[U_truth,S_truth,V_truth] = svd(X_obj);



%% EXAMINE COMPONENTS
if 1
    O = subplot(1,1,1);
    P = plot(diag(S_truth)./S_truth(1,1));
    set(P,'LineWidth',3)
    set(O,'FontSize',20)
    T=title('Ground Truth Scaled Singular Values')
    set(T,'FontSize',20)
    axis([0 25 0 1])
    print -dpng pl_scaledSingVal_GroundTruth

    plot(V_truth(:,1:4))
    title('First 4 times series from Ground Truth SVD')
    print -dpng pl_timeVecSVD_GroundTruth
    
    for jj = 1:4
        subplot(2,2,jj)
        imagesc(reshape(U_truth(:,jj),[N N]))
        colormap gray
    end
    title('First 4 spatial maps from Ground Truth SVD')
    print -dpng pl_spatMapSVD_GroundTruth
    
end



%% DATA SIMULATION

% Set up k-space sampling for navigator and imaging data 
kx_Nav = repmat([-N/2:N/2-1],[1 N_tpt]);    % First navigator is just a line through k-space
ky_Nav = zeros(N,N_tpt,1);
kx_Nav = reshape(kx_Nav,[N N_tpt]);
ky_Nav = reshape(ky_Nav,[N N_tpt]);

% Imaging gradients will read out one line of k-space for every
% acquisition, but will randomize the ky order across each repetition
% FIRST JUST LINEAR ORDERING
kx_image = repmat([-N/2:N/2-1],[1 N N_frames]);   % we will read out one line in k-space for every acquisition
ky_image = repmat([-N/2:N/2-1],[N 1 N_frames]);
kx_image = reshape(kx_image,[N N_tpt]);
ky_image = reshape(ky_image,[N N_tpt]);


if 1
    % INSTEAD, we will use random phase encode ordering
kx_image = repmat([-N/2:N/2-1],[1 N N_frames]);   % we will read out one line in k-space for every acquisition
ky_image = [];
for ii = 1:N_frames
    ky_image = [ky_image,randperm(N)'];
end
ky_image = ky_image - N/2 -1;
ky_image = repmat(ky_image(:)',[N 1]);
kx_image = reshape(kx_image,[N N_tpt]);
ky_image = reshape(ky_image,[N N_tpt]);
end

%% PLOT of K-SPACE TRAJECTORY OVER TIME
J = subplot(1,1,1);
K = plot(ky_image(1,:),'r.');
set(K,'MarkerSize',12);
hold on
M = plot(ky_Nav(1,:),'b.');
set(M,'MarkerSize',12);
set(J,'FontSize',18);
hold off
LL = legend('Imaging','Navigator');
set(LL,'FontSize',18);
print -dpng SamplingPattern

%% SIMULATION OF SAMPLING K-SPACE

% Set up for simulation
[y_im,x_im] = ndgrid(([1:N]-N/2)./N);  % image space goes -1/2 to 1/2 FOV

% Simulation NAVIGATOR SIGNAL == JUST GETTING the CENTER LINE OF KSPACE
% EACH TIME
for ii = 1:N_tpt
    for jj = 1:N
        sig_Nav(jj,ii) = (1/N^2)*sum(sum(obj(:,:,ii).*exp(-j*2*pi*(repmat(kx_Nav(jj,ii),[N N]).*x_im + repmat(ky_Nav(jj,ii),[N N]).*y_im))));% 
    end
end

% Simulation IMAGING SIGNAL == GET RANDOMIZED OR LINEAR LINE OF KSPACE
for ii = 1:N_tpt
    for jj = 1:N
        sig_IMAGING(jj,ii) = (1/N^2)*sum(sum(obj(:,:,ii).*exp(-j*2*pi*(repmat(kx_image(jj,ii),[N N]).*x_im + repmat(ky_image(jj,ii),[N N]).*y_im))));% 
    end
end



%% PERFORM NYQUIST BASED FOURIER TRANSFORM RECONSTRUCTION
% Overlapped sliding window reconstruction
window_wid = 64;
num_imgs = N_tpt - window_wid + 1;
ky_ind = ky_image(1,:)+N/2+1;

for jj = 1:num_imgs
    data_recon=zeros(N,N);
    for ll = 1:window_wid
        data_recon(ky_ind(jj+ll-1),:) = sig_IMAGING(:,jj+ll-1);
    end
    im_Fourier(:,:,jj) = fftshift(ifft2(fftshift(data_recon)));
end


figure
if flag_make_movies    % WRITE OUT MOVIE
H = VideoWriter('Fourier_SlideWind','MPEG-4');
H.FrameRate=10;
open(H)
for ii = 1:num_imgs
    imagesc(abs(im_Fourier(:,:,ii)))
    colormap gray
    %pause(0.1)
    currframe=getframe;
    writeVideo(H,currframe);
end
close(H)
end

%% RECONSTRUCTION WITH LOW RANK APPROACH

% First get time series basis from Navigator data
[U_Nav,S_Nav,V_Nav] = svd(sig_Nav);  % TAKE SVD OF NAVIGATOR DATA
L = 10;  % Rank used
temp_basis = V_Nav(:,1:L);   %% THIS IS THE TEMPORAL BASIS FOR OUR LOW RANK MODEL

% WE DEFINE A LOW RANK STRUCTURE THAT ASSISTS IN OUR ITERATIVE RECONSTRUCTION
A = lr_model(temp_basis,ky_image(1,:)+N/2+1,N); 

spat_map_init = zeros(N, N, L);
niter = 20;
 mask = ones(N*N*N_tpt,1);
 
 % ESTIMATE THE SPATIAL MAPS ASSOCIATED WITH THE TEMPORAL BASES.
 [xs, info] = qpwls_pcg(spat_map_init(:), A, 1, sig_IMAGING(:), 0, 0, 1, niter, mask, 0);
UU = reshape(xs(:,end),[N*N L]); % THESE ARE THE SPATIAL BASES

% NOW COMBINE SPATIAL AND TEMPORAL BASES TO GET FULL MOVIE
movie_im = UU*temp_basis.';
movie_im = abs(reshape(movie_im,[N N N_tpt]));
imagesc(movie_im(:,:,100))
colormap gray

figure
if flag_make_movies    % WRITE OUT MOVIE
H = VideoWriter('Recon_L10','MPEG-4');
H.FrameRate=10;
open(H)
for ii = 1:N_tpt
    imagesc(movie_im(:,:,ii))
    colormap gray
    %pause(0.1)
    currframe=getframe;
    writeVideo(H,currframe);
end
close(H)
end


%% PLOTTING A FEW ITEMS RELATED TO THE LOW RANK RECONSTRUCTION

figure
 O = subplot(1,1,1);
    P = plot(abs(temp_basis(1:200,1:4)));
    set(P,'LineWidth',3);
    set(O,'FontSize',20);
    T=title('First Four Temporal Bases Est. from Navigator Data');
    set(T,'FontSize',20);
    LL = legend('1','2','3','4');
    set(LL,'FontSize',20);
    print -dpng Est_firstFourTempBases
    
    
    figure
    UU = reshape(UU, [N, N, L]);
    subplot(2,2,1),imagesc(abs(UU(:,:,1))); colormap gray
    subplot(2,2,2),imagesc(abs(UU(:,:,2))); colormap gray
    subplot(2,2,3),imagesc(abs(UU(:,:,3))); colormap gray
    subplot(2,2,4),imagesc(abs(UU(:,:,4))); colormap gray
    T = suptitle('First Four Spatial Basis From Imaging Data');
    set(T,'FontSize',20)
    print -dpng Est_firstFourSpatBases
 
