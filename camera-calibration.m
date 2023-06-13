clc;
clear all;
close all;

size_checkerboard = 130; %milimeter
size_images = 5;
number_points = 9;

I = [];
target_points = [];
H = {};
folder_path = 'FixedCamera_Data2/';

for i=1:size_images
    file_name = sprintf('board_%d.jpeg', i);
    file = fullfile(folder_path, file_name);
    img = imread(file);
    img = imresize(img, [3024, 4032]);
    I(:, : ,i) = rgb2gray(img);
end

I = uint8(I);
%%
[coords, ima_pattern] = get_real_points_checkerboard_vmmc(number_points, size_checkerboard, 1);
%%
for i=1:size_images
    image = I(:, :, i);
    target_points(:, :, i) = get_user_points_vmmc(image);
end
%%
transformed = {size(I,3)};
for i=1:size_images
    target_point = target_points(:, :, i);
    H{i} = homography_solve_vmmc(coords', target_point);
    tform = maketform('projective', H{i}');
    transformed{i} = imtransform(ima_pattern, tform);
    figure(i);
    imshow((transformed{i}));
    title('Restored Image');
end

%% Zang's method
H_refined = {};
new_output_img = {};
for i=1:5
    target_point = target_points(:, :, i);
    H_refined{i} = homography_refine_vmmc(coords', target_point, H{i});
    tfrom_new = maketform('projective', H_refined{i}');
    new_output_img{i} = imtransform(ima_pattern, tfrom_new);
    figure(i);
    imshow(new_output_img{i});
    title('New Restored Image');
end

%% get the intiristic parameters
A = internal_parameters_solve_vmmc(H);

%%
save unit_1.mat