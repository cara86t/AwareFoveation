function [sigma_map] = run_on_image(im, sfc, gaze_position)
%RUN_ON_IMAGE Runs the predictor on single image

path_saved = path;
addpath(genpath(fullfile(pwd, 'mex')));
addpath(genpath(fullfile(pwd, 'matlabpyrtools')));

im = im2double(im);
sz = size(im);

% Load display parameters
load(fullfile(pwd, 'disp1_luminance.mat')); % calibration data for pixel intensity to luminance conversion
load(fullfile(pwd, 'disp2_luminance.mat')); % they load disp1_lum and disp2_lum variables
if all(sz == [1440 2560 3])
    % Load the parameters for the 4K display
    disp_prm = display_params(1);
    dlum = disp1_lum;
elseif all(sz == [2160 3840 3])
    % Load the parameters for the QHD display
    disp_prm = display_params(2);
    dlum = disp2_lum;
else
    error('Display parameters are not defined for this input size.');
end
d = disp_prm.distanceToScreen;
s = [disp_prm.screenWidth disp_prm.screenHeight];
r = [disp_prm.resolutionHorizontal disp_prm.resolutionVertical];

% Convert image to luminance
lum = dlum.convert(im); % when using a different display, you will need your own routine for this conversion
fprintf('Luminance range (min, max): %g, %g\n', min(lum(:)), max(lum(:)));

% Load the learned parameters of the predictor for the given subsampling factor
params = get_params(sfc);
alpha = params.alpha;
beta = params.beta;
fec = params.fec;
x1 = params.x1;
x2 = params.x2;
x3 = params.x3;
x4 = params.x4;
omega = params.omega;

% Compute the eccentricity map
szpx = mean(s./r); % pixel size
distpx = d ./ szpx; % distance to the screen in pixels
gaze_coord = gaze_position .* r;
[x,y] = meshgrid(1:r(1), 1:r(2));
x = x - gaze_coord(1);
y = y - gaze_coord(2);
eucdist = sqrt(x.^2 + y.^2);
ecc_map = atand(eucdist ./ distpx);

func = @(block_struct) main(sfc, ...
    block_struct.data(:,:,1), ... % luminance
    mean(mean(block_struct.data(:,:,2))), ... % average eccentricity for the patch
    d, ...
    s, ...
    r, ...
    alpha, ...
    beta, ...
    fec, ...
    x1, ...
    x2, ...
    x3, ...
    x4, ...
    omega);

input = cat(3, lum, ecc_map);
patch_size = [128 128];
overlap_amount = 0.5;
border_size = floor(patch_size .* overlap_amount / 2 );
patch_size = patch_size - 2*border_size;
sigma_map = blockproc(input, ...
    patch_size, ...
    func, ...
    'BorderSize', border_size, ...
    'TrimBorder', false, ...
    'PadMethod', 'symmetric', ...
    'PadPartialBlocks', true, ...
    'UseParallel', false, ...
    'DisplayWaitbar', true);
sigma_map = imresize(sigma_map, [sz(1) sz(2)], 'cubic');

path(path_saved);

end

