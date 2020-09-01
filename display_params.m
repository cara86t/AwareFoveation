function [ r ] = display_params( id )
%DISPLAY_PARAMS
%   Returns the display parameters for the experimental setup

switch id
    case 1
        % Display 1 (QHD)
        r.resolutionVertical = 1440;
        r.resolutionHorizontal = 2560;
        r.distanceToScreen = 66.5;
        r.screenWidth = 59.674;
        r.screenHeight = r.screenWidth / r.resolutionHorizontal * r.resolutionVertical;
        r.blurgamma = 1.8; % applied before Gaussian blur
    case 2
        % Display 2 (4K)
        r.resolutionVertical = 2160;
        r.resolutionHorizontal = 3840;
        r.distanceToScreen = 71;
        r.screenWidth = 69.7829;
        r.screenHeight = r.screenWidth / r.resolutionHorizontal * r.resolutionVertical;
        r.blurgamma = 2.1; % applied before Gaussian blur
end

end

