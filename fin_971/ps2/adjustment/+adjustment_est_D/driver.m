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
M_.fname = 'adjustment_est_D';
M_.dynare_version = '4.8-unstable-2021-11-09-1819-43ac7633';
oo_.dynare_version = '4.8-unstable-2021-11-09-1819-43ac7633';
options_.dynare_version = '4.8-unstable-2021-11-09-1819-43ac7633';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'eps_z'};
M_.exo_names_tex(1) = {'eps\_z'};
M_.exo_names_long(1) = {'eps_z'};
M_.exo_names(2) = {'measure_error'};
M_.exo_names_tex(2) = {'measure\_error'};
M_.exo_names_long(2) = {'measure_error'};
M_.endo_names = cell(8,1);
M_.endo_names_tex = cell(8,1);
M_.endo_names_long = cell(8,1);
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
M_.endo_names(7) = {'qob'};
M_.endo_names_tex(7) = {'qob'};
M_.endo_names_long(7) = {'qob'};
M_.endo_names(8) = {'kob'};
M_.endo_names_tex(8) = {'kob'};
M_.endo_names_long(8) = {'kob'};
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
M_.param_names(4) = {'psi'};
M_.param_names_tex(4) = {'psi'};
M_.param_names_long(4) = {'psi'};
M_.param_names(5) = {'rho'};
M_.param_names_tex(5) = {'rho'};
M_.param_names_long(5) = {'rho'};
M_.param_names(6) = {'gamma'};
M_.param_names_tex(6) = {'gamma'};
M_.param_names_long(6) = {'gamma'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 8;
M_.param_nbr = 6;
M_.orig_endo_nbr = 8;
M_.aux_vars = [];
options_.varobs = cell(2, 1);
options_.varobs(1)  = {'k'};
options_.varobs(2)  = {'qob'};
options_.varobs_id = [ 1 7  ];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
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
M_.nonzero_hessian_eqs = [1 2 5 6];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 8;
M_.eq_nbr = 8;
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
 0 5 12;
 0 6 13;
 2 7 14;
 0 8 15;
 3 9 0;
 0 10 0;
 0 11 0;]';
M_.nstatic = 2;
M_.nfwrd   = 3;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 4;
M_.nspred   = 3;
M_.ndynamic   = 6;
M_.dynamic_tmp_nbr = [7; 5; 1; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'q' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'k' ;
  4 , 'name' , 'z' ;
  5 , 'name' , 'sdf' ;
  6 , 'name' , 'c' ;
  7 , 'name' , 'qob' ;
  8 , 'name' , 'kob' ;
};
M_.mapping.k.eqidx = [1 2 3 6 8 ];
M_.mapping.q.eqidx = [1 2 7 ];
M_.mapping.i.eqidx = [1 2 3 6 ];
M_.mapping.z.eqidx = [2 4 6 ];
M_.mapping.sdf.eqidx = [2 5 ];
M_.mapping.c.eqidx = [5 6 ];
M_.mapping.qob.eqidx = [7 ];
M_.mapping.kob.eqidx = [8 ];
M_.mapping.eps_z.eqidx = [4 ];
M_.mapping.measure_error.eqidx = [7 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 4 6 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(8, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(6, 1);
M_.endo_trends = struct('deflator', cell(8, 1), 'log_deflator', cell(8, 1), 'growth_factor', cell(8, 1), 'log_growth_factor', cell(8, 1));
M_.NNZDerivatives = [27; 27; -1; ];
M_.static_tmp_nbr = [6; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(4) = 1;
oo_.steady_state(2) = 1;
oo_.steady_state(1) = 70;
oo_.steady_state(3) = 10;
oo_.steady_state(5) = 0.96;
oo_.steady_state(6) = 9;
oo_.steady_state(8) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 2) = 0*sqrt(M_.Sigma_e(1, 1)*M_.Sigma_e(2, 2));
M_.Sigma_e(2, 1) = M_.Sigma_e(1, 2);
M_.Correlation_matrix(1, 2) = 0;
M_.Correlation_matrix(2, 1) = M_.Correlation_matrix(1, 2);
M_.sigma_e_is_diagonal = 0;
M_.params(2) = 0.04;
r = M_.params(2);
M_.params(3) = 0.15;
delta = M_.params(3);
M_.params(5) = 0.7;
rho = M_.params(5);
M_.params(6) = 2.0;
gamma = M_.params(6);
estim_params_.var_exo = zeros(0, 10);
estim_params_.var_endo = zeros(0, 10);
estim_params_.corrx = zeros(0, 11);
estim_params_.corrn = zeros(0, 11);
estim_params_.param_vals = zeros(0, 10);
estim_params_.param_vals = [estim_params_.param_vals; 1, NaN, (-Inf), Inf, 5, NaN, NaN, 0, 1, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, NaN, (-Inf), Inf, 5, NaN, NaN, 0, 1, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 1, NaN, (-Inf), Inf, 4, 78.1330, Inf, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 2, NaN, (-Inf), Inf, 4, 0.1624, Inf, NaN, NaN, NaN ];
options_.mh_jscale = 1.8;
options_.mode_compute = 6;
options_.datafile = 'data_ps2';
options_.order = 1;
var_list_ = {'k';'qob'};
oo_recursive_=dynare_estimation(var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'adjustment_est_D_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'adjustment_est_D_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'adjustment_est_D_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'adjustment_est_D_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'adjustment_est_D_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'adjustment_est_D_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'adjustment_est_D_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
