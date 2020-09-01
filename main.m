function rs = main( sfc, im, vd, d, s, r, alpha, beta, fec, x1, x2, x3, x4, omega )
%MAIN
% Estimates the tolereable amount of blur in the peripheral visual field
% (as Stdev of Gaussian kernel) for the eccentricity specified by input
% argument vd.
%
%   INPUTS 
%
%   sfc: Subsampling factor
%
%   patch: Reference patch (single channel image with absolute luminance in
%   cd/m^2)
%
%   vd: Eccentricity (visual degrees)
%
%   d: The distance of the observer to the screen (scalar in cm)
%
%   s: An array of length 2 which contains the horizontal and vertical size
%   of the screen in centimeters
%
%   r: An array of length 2 which contains the horizontal and vertical
%   screen resolution
%
%   alpha: Self-masking parameter
%
%   beta:  Spatial-masking parameter
%
%   fec:   Fundamental eccentricity constant. It controls how fast the CSF
%   declines with eccentricity
%
%   x1, x2, x3, x4: CSF parameters
%
%   omega: Softmax parameter for pooling
%   
%	OUTPUTS
%
%	rs:  The estimated Stdev of Gaussian kernel

if sfc > 1
    im = im(1:sfc:end, 1:sfc:end); % comment this line out if the input is already subsampled
end

% pixel size
psz = s ./ r;
pszm = mean(psz); % inter-pixel distance assuming square pixels

% peak CPD supported by the display for the given eccentricity
pcpd = d .* (tand(vd+0.5) - tand(vd-0.5)) ./ pszm ./ 2 ./ sfc;

decomp_type = 'laplacian'; % use only laplacian for now
a = 0.4;
filt = [1/4-a/2 1/4 a 1/4 1/4-a/2].'; % filter from Burt's [1983] paper

switch decomp_type
    case 'laplacian'
        % pyr_height = max(min(floor(log2(pcpd))+2, 6), 3); % pyrtools supports maximum 7 levels at 2k res
        switch sfc
        case 1
            pyr_height = 6;
        case 2
            pyr_height = 5;
        case 4
            pyr_height = 4;
        otherwise
            error('pyr_height is not defined for the given subsampling factor (sfc).');
        end
        filt_l = sqrt(2) .* filt;
        [lpyr, lpind] = buildLpyr(im, pyr_height, filt_l);
    case 'steerable'
        % this is not tested yet (do not use)
        pyr_height = min(floor(log2(pcpd))+3, 5); % pyrtools supports maximum 5 levels
        filts = 'sp5Filters';
        [lpyr,lpind] = buildSpyr(im, pyr_height, filts);  
end

peak_pix_freq = 2.^(-(1:pyr_height)) / sfc; % peak frequencies in terms of pixels
peak_freq = pcpd .* 2.^(-(1:pyr_height)) * 2;

[gpyr, gpind] = buildGpyr(im, pyr_height, filt);

% Compute the band-limited Weber contrasts and the contribution of multiple
% bandpass mechanisms
weber_contrast = cell([pyr_height-2 1]);
contrast_sign = cell([pyr_height-2 1]);
local_luminance = cell([pyr_height-2 1]);
jnd = cell([pyr_height-2 1]);
pass_coeff = cell([pyr_height-2 1]);
F_csf = cell([pyr_height 1]);
self_masking = cell([pyr_height-2 1]);
spatial_masking = cell([pyr_height-2 1]);
mask_filter = ones(5,5);
mask_filter(3,3) = 0;
normalization_c = 1./(sum(mask_filter(:)));
mask_filter = normalization_c .* mask_filter;

for i = 1:pyr_height-2
    lum_diff = pyrBand(lpyr, lpind, i);

    local_luminance{i} = pyrBand(gpyr, gpind, i+2);
    local_luminance{i} = upConv(local_luminance{i}, filt * filt.' * 4, 'reflect1', [2 2]);
    local_luminance{i} = upConv(local_luminance{i}, filt * filt.' * 4, 'reflect1', [2 2]);
    local_luminance{i} = local_luminance{i}(1:size(lum_diff,1),1:size(lum_diff,2));
    weber_contrast{i} = lum_diff ./ local_luminance{i};
    weber_contrast{i}(isnan(weber_contrast{i})) = 0; % NaNs may result from 0/0
    contrast_sign{i} = sign(weber_contrast{i});

    % Cubic interpolation between the sensitivities of 4 bands at 4, 8, 16 and 32 cpds
    log_s = interp1(log([1e-32 4 8 16 32 64]), ...
        log(horzcat(10.^[x1 x1 x2 x3 x4], [1e-32])), ...
        log(peak_freq(i)), 'pchip', -Inf);

    barten = @(lum) power(1 + 0.7 ./ lum, -0.2);
    F_csf{i} = exp(log_s) ./ exp(fec.*vd.*peak_freq(i)) .* barten(local_luminance{i});

    jnd{i} = abs(weber_contrast{i} .* F_csf{i});
    
    self_masking{i} = power(jnd{i}, alpha);
    spatial_masking{i} = imfilter(power(abs(jnd{i}), beta), mask_filter, 'symmetric');
    tmp = self_masking{i} - 1 - spatial_masking{i};
    pass_coeff{i} = power(tmp, 1./alpha) ./ jnd{i};
    pass_coeff{i}(tmp < 0) = 0; % if the masked signal is below 1 JND
end

pass_coeff_mat = zeros([size(im) pyr_height-2]);
for i = 1:numel(pass_coeff)
    pass_coeff_mat(:,:,i) = imresize(pass_coeff{i}, size(im), ...
        'Method', 'bilinear', ...
        'Antialiasing', false, ...
        'Dither', false);
end
sz = size(pass_coeff_mat);
% Compute sigma estimate for each pixel
if numel(sz) == 2 % pass coeff mat has only one level, this occurs if the input patch is very small or the subsampling factor is large
    sz = [sz 1];
end

frequency_bands = repmat(reshape(peak_pix_freq(1:pyr_height-2), [1 1 sz(3)]), sz(1:2));
sigma_f = sqrt(-power(frequency_bands,2) ./ (2 .* log(pass_coeff_mat)));
sigma_s = 1 ./ (2 * pi * sigma_f);
[sigma_est, ~] = min(sigma_s,[],3); % sigma_est contains the estimations for each pixel (before pooling)

% Pooling
pooling_param = sign(omega) * (2.^(omega)-1);
rl = 1./sigma_est;

% Apply smooth max in spatial domain
rs_exp = exp(pooling_param .* rl);
rs_nom = rl .* rs_exp;
rs_nom = sum(rs_nom(:));
if isinf(rs_nom)
    rs = max(rl(:));
elseif rs_nom == 0
    rs = min(rl(:));
else
    rs_denom = sum(rs_exp(:));
    rs = rs_nom ./ rs_denom;
end
rs = 1 ./ rs;
rs = min(rs, 10);



end
