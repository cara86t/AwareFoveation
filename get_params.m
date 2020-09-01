function [ params ] = get_params( sfc )
%GET_PARAMS Returns optimal parameter values learned from Simulated
%Annealing

filename = fullfile(pwd, ...
    'matfiles', ...
    sprintf('best_params_%d.mat', sfc));
if ~exist(filename, 'file')
    error('Get_params: cannot read the file: %s', filename);
%     params.alpha = 0.5;
%     params.beta = 0.2;
%     params.fec = 0.035;
%     params.omega = 4;
else
    p = load(filename);
    params = p.params;
    fprintf('Loaded parameters from %s\n', ...
        filename);
end


end

