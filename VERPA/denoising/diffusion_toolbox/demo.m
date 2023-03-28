% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

clear
close all
run('set_path.m');
img = imread('img/monet_bridge.png');

%% computing the filtered images
options = get_diff_options;
f1_pm = perona_malik_diff(img, options);

options = get_aniso_diff_options();
options.diffusivity = 'weickert';
f1_edge = aniso_diff_filt(img, options);

options = get_aniso_coh_diff_options();
f1_coh = aniso_coh_diff_filt(img, options);

%% plotting:
subplot(2, 4, 1);
imshow(img);
title('raw frame');

subplot(2, 4, 2);
imshow(f1_pm);
title('nonlinear, isotropic diffusion');

subplot(2, 4, 3);
imshow(f1_edge);
title('edge enhancing diffusion');

figure(1);
subplot(2, 4, 4);
imshow(f1_coh);
title('coherence enhancing diffusion');

bb = 350;
subplot(2, 4, 5);
imshow(img(1:bb, end-bb-1:end, :));
subplot(2, 4, 6);
imshow(f1_pm(1:bb, end-bb-1:end, :));
subplot(2, 4, 7);
imshow(f1_edge(1:bb, end-bb-1:end, :));
figure(1);
subplot(2, 4, 8);
imshow(f1_coh(1:bb, end-bb-1:end, :));
