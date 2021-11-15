%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'noadjustment';
M_.dynare_version = '4.8-unstable-2021-11-09-1819-43ac7633';
oo_.dynare_version = '4.8-unstable-2021-11-09-1819-43ac7633';
options_.dynare_version = '4.8-unstable-2021-11-09-1819-43ac7633';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps_z'};
M_.exo_names_tex(1) = {'eps\_z'};
M_.exo_names_long(1) = {'eps_z'};
M_.endo_names = cell(6,1);
M_.endo_names_tex = cell(6,1);
M_.endo_names_long = cell(6,1);
M_.endo_names(1) = {'k'};
M_.endo_names_tex(1) = {'k'};
M_.endo_names_long(1) = {'k'};
M_.endo_names(2) = {'q'};
M_.endo_names_tex(2) = {'q'};
M_.endo_names_long(2) = {'q'};
M_.endo_names(3) = {'i'};
M_.endo_names_tex(3) = {'i'};
M_.endo_names_long(3) = {'i'};
M_.endo_names(4) = {'z'};
M_.endo_names_tex(4) = {'z'};
M_.endo_names_long(4) = {'z'};
M_.endo_names(5) = {'sdf'};
M_.endo_names_tex(5) = {'sdf'};
M_.endo_names_long(5) = {'sdf'};
M_.endo_names(6) = {'c'};
M_.endo_names_tex(6) = {'c'};
M_.endo_names_long(6) = {'c'};
M_.endo_partitions = struct();
M_.param_names = cell(6,1);
M_.param_names_tex = cell(6,1);
M_.param_names_long = cell(6,1);
M_.param_names(1) = {'theta'};
M_.param_names_tex(1) = {'theta'};
M_.param_names_long(1) = {'theta'};
M_.param_names(2) = {'r'};
M_.param_names_tex(2) = {'r'};
M_.param_names_long(2) = {'r'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'delta'};
M_.param_names_long(3) = {'delta'};
M_.param_names(4) = {'rho'};
M_.param_names_tex(4) = {'rho'};
M_.param_names_long(4) = {'rho'};
M_.param_names(5) = {'sigma'};
M_.param_names_tex(5) = {'sigma'};
M_.param_names_long(5) = {'sigma'};
M_.param_names(6) = {'gamma'};
M_.param_names_tex(6) = {'gamma'};
M_.param_names_long(6) = {'gamma'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 6;
M_.param_nbr = 6;
M_.orig_endo_nbr = 6;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [2 5 6];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 6;
M_.eq_nbr = 6;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 4 0;
 0 5 0;
 0 6 0;
 2 7 10;
 0 8 11;
 3 9 0;]';
M_.nstatic = 2;
M_.nfwrd   = 1;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 2;
M_.nspred   = 3;
M_.ndynamic   = 4;
M_.dynamic_tmp_nbr = [3; 4; 1; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'q' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'k' ;
  4 , 'name' , 'z' ;
  5 , 'name' , 'sdf' ;
  6 , 'name' , 'c' ;
};
M_.mapping.k.eqidx = [2 3 6 ];
M_.mapping.q.eqidx = [1 ];
M_.mapping.i.eqidx = [3 6 ];
M_.mapping.z.eqidx = [2 4 6 ];
M_.mapping.sdf.eqidx = [2 5 ];
M_.mapping.c.eqidx = [5 6 ];
M_.mapping.eps_z.eqidx = [4 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 4 6 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(6, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(6, 1);
M_.endo_trends = struct('deflator', cell(6, 1), 'log_deflator', cell(6, 1), 'growth_factor', cell(6, 1), 'log_growth_factor', cell(6, 1));
M_.NNZDerivatives = [17; 14; -1; ];
M_.static_tmp_nbr = [2; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 0.7;
theta = M_.params(1);
M_.params(2) = 0.04;
r = M_.params(2);
M_.params(3) = 0.15;
delta = M_.params(3);
M_.params(4) = 0.7;
rho = M_.params(4);
M_.params(5) = 0.01;
sigma = M_.params(5);
M_.params(6) = 2.0;
gamma = M_.params(6);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(4) = 1;
oo_.steady_state(2) = 1;
oo_.steady_state(1) = (M_.params(1)/(M_.params(2)+M_.params(3)))^(1/(1-M_.params(1)));
oo_.steady_state(3) = M_.params(3)*oo_.steady_state(1);
oo_.steady_state(5) = 1/(1+M_.params(2));
oo_.steady_state(6) = oo_.steady_state(1)^M_.params(1)-oo_.steady_state(3);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
resid(1);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(5))^2;
options_.irf = 40;
options_.nodisplay = true;
options_.order = 2;
options_.periods = 200;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'noadjustment_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'noadjustment_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'noadjustment_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'noadjustment_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'noadjustment_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'noadjustment_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'noadjustment_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
