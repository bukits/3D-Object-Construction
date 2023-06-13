clear all
close all
clc
%% PARAMETER DEFINITION %%
%%define dataset path
params.Directory    = fullfile('dataset');

%%Detector parameters
params.detector     =  'KAZE'; %'LoG_ss', 'SURF', 'SIFT', ...
params.nscales      =        10;
params.noctaves     =        3;
params.sigma0       =      1.6; % as we are using Matlab functions this is the minimum value allowed
params.npoints      =      350;
params.th           =    0.001; % alternative to npoints.
imscale = 0.5;

%%Descriptor parameters (see doc extractFeatures for explanation and additional parameters)
params.descriptor   =  'KAZE';  % 'SIFT', 'SURF', 'DSP-SIFT'...
params.desOnDecom   =    false; % describe on scale-space (linear or non-linear) decomposition (if available)
params.Upright      =   false;  % set to true to avoid orientation estimation.
% for DSP-SIFT
params.dsp.ns       =      15;% number of sampled scales
params.dsp.sc_min   =     1/6;% smallest scale (relative to detection)
params.dsp.sc_max   =       7;% largest scale (relative to detection);    

%%Matching parameters (see doc matchFeatures for explanation and additional parameters)
params.MaxRatio     =   0.5;
params.Metric       =  'SSD';
%% END OF PARAMETER DEFINITION %%
%% addpaths
addpath(genpath('./detectors/'));
addpath(genpath('./descriptors/'));
addpath(genpath('./toolbox/'));
%% preload dataset
params.Scene = imageDatastore(params.Directory);
numImages    = numel(params.Scene.Files);
%% initialize (sort of)
ima{numImages}           = [];
points{numImages}        = [];
decomposition{numImages} = [];
features{numImages}      = [];
%% get sigmas
k = 1;
params.sigmas = zeros(1,params.noctaves*params.nscales);
for o = 0:params.noctaves-1
    params.sigmas(k:(k+params.nscales-1)) = params.sigma0.*pow2([0:(params.nscales-1)]/params.nscales + o);
    k = k+params.nscales;
end
%% detect & describe
for j = 1:numImages
    %% Load and convert images %%
    ima{j} = readimage(params.Scene, j);
    ima{j} = imresize(ima{j}, imscale); % resizing the image 0.5 of the original size
    gima = im2double(rgb2gray(ima{j}));
    %% PoI Detection %%
    sprintf('Detecting for image: %d',j)
    [points{j},decomposition{j}] = myDetector(gima,params);
    %% PoI Description %%
    sprintf('Describing for image: %d',j)
    [features{j},points{j}] = myDescriptor(points{j},decomposition{j},params);
    %% show detections
    figure(j)
    imshow(ima{j}); hold on;
    plot(points{j},'showOrientation',true);
end
%% montage
figure();
montage(ima);
%% PoI Matching (assumes two images, i.e. numImages == 2) %%
%selecting Image pair
i1=8; i2=9;

indexPairs       = matchFeatures(features{i1},features{i2},'MaxRatio',params.MaxRatio,'Metric',params.Metric) ;
matchedPoints{1} = points{i1}(indexPairs(:,1));
matchedPoints{2} = points{i2}(indexPairs(:,2));
figure(numImages+1); showMatchedFeatures(ima{i1},ima{i2},matchedPoints{1},matchedPoints{2});
legend('matched points 1','matched points 2');
%% Homography estimation and warp %% 
%% A) Estimate the transformation between ima(2) and ima(1).
if numel(matchedPoints{2}.Scale) < 4
   sprintf('Unable to match enough points -> End of program')
   return;
end
[tform21,in1,in2]  = estimateGeometricTransform(matchedPoints{2}, matchedPoints{1},...
         'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
        
warpedImage = imwarp(ima{i2}, tform21, 'OutputView', imref2d(size(ima{i1})));
% show results
ima2         = zeros(size(ima{i1}));
for ch=1:3
ima2(:,:,ch) = imresize(ima{i2}(:,:,ch),size(ima{i1}(:,:,ch)));
end
multi = cat(4,ima{i1},ima2,ima{i1},warpedImage);
figure(numImages+2);aa = montage(multi,'Size',[2,2]);
result21 = aa.CData;
disp  = 20;
figure(numImages+2);clf,imshow(result21)
text(disp,disp,'Image 1','Color','red','FontSize',14)
text(disp + size(result21,2)/2,disp,'Image 2','Color','red','FontSize',14)
text(disp,disp + size(result21,1)/2,'Image 1','Color','red','FontSize',14)
text(disp + size(result21,2)/2,disp + size(result21,1)/2,'Image 2 to 1','Color','red','FontSize',14)
%% B) Estimate the transformation between ima(1) and ima(2).
[tform12, in3, in4]  = estimateGeometricTransform(matchedPoints{1}, matchedPoints{2},...
         'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
        
warpedImage = imwarp(ima{i1}, tform12, 'OutputView', imref2d(size(ima{i2})));
% show results
for ch=1:3
ima1(:,:,ch) = imresize(ima{i1}(:,:,ch),size(ima{i2}(:,:,ch)));
end
multi = cat(4,ima1,ima{i2},warpedImage,ima{i2});
figure(numImages+3);aa = montage(multi,'Size',[2,2]);
result12 = aa.CData;
figure(numImages+3);clf,imshow(result12)
text(disp,disp,'Image 1','Color','red','FontSize',14)
text(disp + size(result12,2)/2,disp,'Image 2','Color','red','FontSize',14)
text(disp,disp + size(result12,1)/2,'Image 1 to 2','Color','red','FontSize',14)
text(disp + size(result12,2)/2,disp + size(result12,1)/2,'Image 2','Color','red','FontSize',14)
%% Quantitative and qualitative evaluationa
[F, inliersIndex] = estimateFundamentalMatrix(matchedPoints{1}, matchedPoints{2});

rnk = rank(F);

fprintf('Number of inliers, calculating the homography transformation between the two views is %d points \n', in1.Count);
fprintf('Number of inliers, calculating the fundamental matrix between the two views is %d points \n', sum(inliersIndex));
fprintf('Total number of matched points is %d\n', numel(matchedPoints{2}.Scale));

vgg_gui_F(ima{i2}, ima{i1}, F');

save unit_2.mat


