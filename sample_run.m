function sample_run()
%SAMPLE_RUN Runs the predictor on a sample scene

im_4K = imread(fullfile(pwd, 'scenes', '000.jpg')); % load the original scene (4K)
im_QHD = imresize(im_4K, [1440 2560]); % resize to simulate low-res display

% Setting sfc = 1 here runs the predictor on full-resolution inputs with no
% downsampling. That provides a better fit to human perception but runs much
% slower and it's not very useful in practice because it does not provide any
% speedup in rendering.
% You may consider using sfc = 1 if prediction performance has more 
% priority than the computational cost incurred by having a full-resolution
% input or additional computational costs are not critical for your
% application (e.g. using the predictions as an input in an offline
% processing pipeline).
% sfc = 4 gives similar predictions to those which are used in our subjective
% experiments.
sfc = 1; % subsampling factor
gaze_position = [0.5 0.5]; % assumed to be at the center of the screen
% [0,0] and [1,1] corresponds to the corners of the screen

%% This should reproduce the result for scene G1 in Fig. 10
sigma_map_QHD = run_on_image(im_QHD, sfc, gaze_position);

fig1 = figure('Name', 'Input Image QHD');
imshow(im_QHD);
title(fig1.CurrentAxes, 'Input Image QHD');

fig2 = figure('Name', 'Sigma Map QHD');
imshow(sigma_map_QHD, 'DisplayRange', [0 2.5], 'Colormap', flipud(parula(1000)));
colorbar(fig2.CurrentAxes);
title(fig2.CurrentAxes, '\sigma map QHD');

%% Running the metric for the 4K display
% This should result in higher sigma predictions because 4K display
% covers a larger visual field
sigma_map_4K = run_on_image(im_4K, sfc, gaze_position);

fig1 = figure('Name', 'Input Image 4K');
imshow(im_4K);
title(fig1.CurrentAxes, 'Input Image 4K');

fig2 = figure('Name', 'Sigma Map 4K');
imshow(sigma_map_4K, 'DisplayRange', [0 4], 'Colormap', flipud(parula(1000)));
colorbar(fig2.CurrentAxes);
title(fig2.CurrentAxes, '\sigma map 4K');

end

