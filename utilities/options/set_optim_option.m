function [noptim] = set_optim_option(optim,objfun,list_of_variables,enr)
% This function checks the optim structure and sets default values if needed.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

noptim = optim;

noptim = setoptions(noptim, 'method', 'gradDesc', {'gradDesc','gradDesc_sto','bfgs','sfista'});
noptim = setoptions(noptim, noptim.method, []);


switch lower(noptim.method)
    case 'graddesc'
        noptim.gradDesc = set_gradDesc_defaults(noptim.gradDesc, objfun, list_of_variables, enr);
    case 'graddesc_sto'
        noptim.gradDesc_sto = set_gradDesc_defaults(noptim.gradDesc_sto, objfun, list_of_variables, enr);
    case 'bfgs'
        noptim = setoptions(noptim,'kernel_size_geom_reg', 0);
        noptim.bfgs = set_bfgs_defaults(noptim.bfgs, list_of_variables, enr);
    case 'sfista'
        noptim = setoptions(noptim,'kernel_size_geom_reg', 0);
        noptim.sfista = set_sfista_defaults(noptim.sfista, list_of_variables, enr);       
        
end

end

function opt = set_gradDesc_defaults(opt, objfun, list_of_variables, enr_name)

opt = setoptions(opt, 'step_increase', 1.2);
opt = setoptions(opt, 'step_decrease', .5);
opt = setoptions(opt, 'kernel_size_signal_reg', 0);
opt = setoptions(opt, 'kernel_size_geom_reg', 0);

switch objfun.distance
    case 'kernel'
        opt = setoptions(opt, 'max_nb_iter', 50 * ones(1, size(objfun.kernel_distance.kernel_size_geom,2)));
        nb_run = length(opt.max_nb_iter);
        
        if ~isequal(length(objfun.kernel_distance.kernel_size_geom),nb_run)...
                || ( (strcmp(objfun.kernel_distance.kernel_grass,'gaussian_oriented') || strcmp(objfun.kernel_distance.kernel_grass,'gaussian_unoriented')) && ~isequal(length(objfun.kernel_distance.kernel_size_grass),nb_run) )
            error('All kernel sizes objfun.sigmaXX and optim.max_nb_iter must have the same length');
        end
    
    case 'wasserstein'
        opt = setoptions(opt, 'max_nb_iter', 50 * ones(1, size(objfun.wasserstein_distance.epsilon,2)));

end

opt = setoptions(opt, 'save_template_evolution', 0);
opt = setoptions(opt, 'min_step_size', 1e-10);
opt = setoptions(opt, 'min_fun_decrease', 1e-4);

for i = 1:length(list_of_variables)
    opt = setoptions(opt, ['step_size_',list_of_variables{i}], 'auto');
end

opt = setoptions(opt, 'list_of_variables', list_of_variables);
opt = setoptions(opt, 'enr_name', enr_name);

end

function opt = set_bfgs_defaults(opt, list_of_variables, enr_name)

opt = setoptions(opt, 'nvec', 20); % BFGS memory
opt = setoptions(opt, 'maxit', 50);
opt = setoptions(opt, 'prtlevel', 0);
opt = setoptions(opt, 'normtol', eps);
opt = setoptions(opt, 'tol', eps);
opt = setoptions(opt, 'list_of_variables', list_of_variables);
opt = setoptions(opt, 'enr_name', enr_name);
opt = setoptions(opt, 'record_history', 0); % record the history of variables...

end

function opt = set_sfista_defaults(opt, list_of_variables, enr_name)

opt = setoptions(opt, 'nvec', 20); % BFGS memory
opt = setoptions(opt, 'maxit', 50);
opt = setoptions(opt, 'prtlevel', 0);
opt = setoptions(opt, 'normtol', eps);
opt = setoptions(opt, 'tol', eps);
opt = setoptions(opt, 'list_of_variables', list_of_variables);
opt = setoptions(opt, 'enr_name', enr_name);
opt = setoptions(opt, 'record_history', 0); % record the history of variables...

end


