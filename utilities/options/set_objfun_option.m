function [nobjfun] = set_objfun_option(objfun,init,target)
% This function checks the objfun structure and sets default values if needed.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

nobjfun = objfun;

nobjfun = setoptions(nobjfun, 'distance', 'kernel', {'kernel','wasserstein','pt2pt'});

switch lower(nobjfun.distance)
    case 'kernel'
        nobjfun.kernel_distance = set_kernel_distance_defaults(nobjfun.kernel_distance);
        
    case 'wasserstein'
        nobjfun.wasserstein_distance = set_wasserstein_distance_defaults(nobjfun.wasserstein_distance);
        
    case 'pt2pt'
        if ~isequal(size(target.x,1),size(init,1))
            error('Check source and target sizes to use pt2pt distance');
        end
end

nobjfun = setoptions(nobjfun, 'weight_coef_dist', 40);
nobjfun = setoptions(nobjfun, 'signal_type','vertex');
nobjfun = setoptions(nobjfun, 'data_signal_type','vertex');

end


function obj = set_wasserstein_distance_defaults(obj)

obj = setoptions(obj, 'epsilon');
obj = setoptions(obj, 'niter', 1000);
obj = setoptions(obj, 'tau', 0); % basic sinkhorn, no extrapolation
obj = setoptions(obj, 'rho', Inf); % balanced case
obj = setoptions(obj, 'weight_cost_varifold', [1,.01]); % weight on spatial and orientation distances
obj = setoptions(obj, 'method', 'matlab', {'cuda','matlab'});

end

function obj = set_kernel_distance_defaults(obj)

obj = setoptions(obj, 'distance', 'empty');

if ~strcmp(obj.distance, 'empty')
    switch lower(obj.distance)
        case 'cur'
            obj.kernel_geom = 'gaussian';
            obj.kernel_signal = 'gaussian';
            obj.kernel_grass = 'linear';
        case 'var'
            obj.kernel_geom = 'gaussian';
            obj.kernel_signal = 'gaussian';
            obj.kernel_grass = 'binet';
        case 'varexpo'
            obj.kernel_geom = 'gaussian';
            obj.kernel_signal = 'gaussian';
            obj.kernel_grass = 'gaussian_unoriented';
        otherwise
            obj.distance = 'var';
            obj.kernel_geom = 'gaussian';
            obj.kernel_signal = 'gaussian';
            obj.kernel_grass = 'binet';
            warning('distance : Possible distance are current (''cur'') or varifold (''var'' or ''varexpo''). objfun.distance is set to ''var''.')
    end
       
elseif ~isfield(obj,'kernel_geom') || ~isfield(obj,'kernel_signal') || ~isfield(obj,'kernel_grass')
    error('Please provide a field distance (cur, var or varexpo) or kernel type for geometry, signal and grassmanian');
end

obj = setoptions(obj, 'kernel_size_geom');  % size of the geometric kernel in the data attachment term
obj = setoptions(obj, 'kernel_size_signal');  % size of the functional kernel in the data attachment term
if ~strcmpi(obj.kernel_grass, 'binet') && ~strcmpi(obj.kernel_grass, 'linear')
    obj = setoptions(obj, 'kernel_size_grass');
end

obj = setoptions(obj, 'method', 'matlab', {'cuda','matlab'});

end
