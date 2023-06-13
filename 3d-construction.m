%% loading the data
clc;
clear all;
close all;

ACT_path = './ACT_lite/';
addpath(genpath(ACT_path));

extra_funs_path = './extra_funs/';
addpath(genpath(extra_funs_path));
allfns = './allfns';
addpath(genpath(allfns));

curr_dir='../unit_3';
dir = '../unit_2';
chdir(dir);
unit2 = load("unit_2.mat");
features = unit2.features;
ima = unit2.ima;
points = unit2.points;
MaxRatio     =      0.8;
Metric       =  unit2.params.Metric;
chdir(curr_dir);
%% import unit 1
dir = '../unit_1';
chdir(dir);
unit1 = load('unit_1.mat');
K = unit1.A;
chdir(curr_dir);
%% create the 3D points
%select cameras
start_range = 4;
end_range = 9;

%1:3 good
images = (start_range:end_range);
ncam = length(images);

features_new = {};
points_new = {};
ima_new = {};

for i=1:ncam
    features_new{i} = features{images(i)};
    points_new{i} = points{images(i)};
    ima_new{i} = ima{images(i)};
end

%% 1 Compute consistent point matches among N views.
q_data = n_view_matching(points_new, features_new, ima_new, MaxRatio, Metric);
q_data = homogenize_coords(q_data);
%%
imscale = unit2.imscale;

K_est = K.*imscale; 
K_est(3, 3) = 1;

K_ = zeros(3, 3, 2);
K_(:, :, 1) = K_est;
K_(:, :, 2) = K_est;
%%
% 2. Compute the Fundamental matrix and an initial projective reconstruction from 2 of the cameras.
% ------------------------------------------------------------------------
q_2cams(:,:,1) = q_data(:,:,1); 
q_2cams(:,:,2) = q_data(:,:,end);
[F, P_2cam_est,Q_2cam_est,q_2cam_est] = MatFunProjectiveCalib(q_2cams);

disp(['Residual reprojection error. 8 point algorithm   = ' num2str( ErrorRetroproy(q_2cams,P_2cam_est,Q_2cam_est)/2 )]);
draw_reproj_error(q_2cams,P_2cam_est,Q_2cam_est);

%%
% ------------------------------------------------------------------------
% 3.a Resection. Obtain the projection matrices of the rest of cameras using the PDLT_NA function 
% ------------------------------------------------------------------------
% ...

P_cams(:,:,:) = zeros(3,4,ncam);
P_cams(:,:,1) = P_2cam_est(:,:,1);
P_cams(:,:,ncam) = P_2cam_est(:,:,2);

for i = 2:ncam-1
    P_cams(:,:,i) = PDLT_NA(q_data(:,:,i),Q_2cam_est);
end

disp(['Residual reprojection error, After resectioning  = ' num2str( ErrorRetroproy(q_data,P_cams,Q_2cam_est)/2 )]);
draw_reproj_error(q_data,P_cams,Q_2cam_est); %q_data

%%
% ------------------------------------------------------------------------
% 3.b Projective Bundle Adjustment. Use BAProjectiveCalib function
% Coordinates of 3D and 2D points are given in homogeneus coordinates
% ------------------------------------------------------------------------
% auxiliary matrix that indicates that all points are visible in all the cameras
npoints = size(q_data,2);
vp = ones(npoints,ncam);
% ...
[P_bundle,Q_bundle,q_bundle] = BAProjectiveCalib(q_data,P_cams,Q_2cam_est,vp);

% ------------------------------------------------------------------------
% Compute the statistics of the reprojection error for the improved projective reconstruction
% ------------------------------------------------------------------------
% ...
disp(['Reprojection error, After Bundle Adjustment  = ' num2str( ErrorRetroproy(q_data,P_bundle,Q_bundle)/2 )]); %q_data
draw_reproj_error(q_data,P_bundle,Q_bundle);

%%
% 4. Obtain the essential matrix (E) from the fundamental matrix (F) and the
% intrinsic parameter matrices (K).
% ------------------------------------------------------------------------

% estimating F from P_BA (P_bundle)
F_bundle = vgg_F_from_P({P_bundle(:,:,1), P_bundle(:,:,end)});

%% 5.
K_1 = K_(:, :, 1);
K_2 = K_(:, :, 2);
E = K_2' * F_bundle * K_1;

% calculate the T & R matrices
[R_est,T_est] = factorize_E(E);

q_bundle_2cam(:,:,1) = q_bundle(:,:,1);
q_bundle_2cam(:,:,2) = q_bundle(:,:,end);

% ------------------------------------------------------------------------
% Save the 4 solutions (R,t) in the structures Rcam(3,3,cam,sol), T(3,cam,sol),
% where cam indicates the camera number and sol indicates the solution number (1, 2, 3 or 4).
% ------------------------------------------------------------------------
Rcam = zeros(3,3,2,4);
Tcam = zeros(3,2,4);
I = eye(3);

Rcam(:, :, 1, 1) = I;
Rcam(:, :, 1, 2) = I;
Rcam(:, :, 1, 3) = I;
Rcam(:, :, 1, 4) = I;

Tcam(:, 1, :) = 0;

Rcam(:, :, 2, 1) = R_est(:, :, 1);
Rcam(:, :, 2, 2) = R_est(:, :, 1);
Rcam(:, :, 2, 3) = R_est(:, :, 2);
Rcam(:, :, 2, 4) = R_est(:, :, 2);

Tcam(:, 2, 1) = T_est;
Tcam(:, 2, 2) = -T_est;
Tcam(:, 2, 3) = T_est;
Tcam(:, 2, 4) = -T_est;
%% test fundemental matrix
vgg_gui_F(ima_new{1}, ima_new{end}, F_bundle');
%%
% ------------------------------------------------------------------------
% 6. For each solution we obtain an Euclidean solution and we visualize it.
% ------------------------------------------------------------------------
npoints = size(q_data,2);
Q_euc = zeros(4,npoints,2); % Variable for recontructed points
P_euc = zeros(3,4,2);       % Variable for projection matrices
figNo=figure;
for sol=1:4
    % Euclidean triangulation to obtain the 3D points (use TriangEuc)
    Q_euc = TriangEuc(Rcam(:, :, 2, sol), Tcam(:, 2, sol), K_, q_bundle_2cam);
       
    % visualize 3D reconstruction
    figure();
    draw_scene(Q_euc, K_, Rcam(:,:,:,sol), Tcam(:,:,sol));
    title(sprintf('Solution %d', sol));
     
    % Compute the projection matrices from K, Rcam, Tcam
    for k=1:2
        P_euc(:, :, k) = K_(:, :, k) * [Rcam(:, :, k, sol) (-1) * Rcam(:, :, k, sol) * Tcam(:, k, sol)];
    end
    
    % Obtain the re-projected points q_rep
    for k=1:2
        q_rep(:, :, k) = P_euc(:, :, k) * Q_euc;
    end
    
    % Visualize reprojectd points to check that all solutions correspond to
    % the projected images
    q_rep = un_homogenize_coords(q_rep);
    for k=1:2
      figure(figNo); subplot(4,2,2*(sol-1)+k); scatter(q_rep(1,:,k),q_rep(2,:,k),30,[1,0,0]);
      title(sprintf('Reprojection %d, image %d', sol, k));
      daspect([1, 1, 1]);
      pbaspect([1, 1, 1]);
      axis([-1000, 1000, -1000, 1000]);
    end
end

%%
disp(['Reprojection error, Euclidean Re-construction  = ' num2str( ErrorRetroproy(q_2cams,P_euc,Q_euc)/2 )]); %q_data
draw_reproj_error(q_2cams,P_euc,Q_euc);

disp('************************************* END')