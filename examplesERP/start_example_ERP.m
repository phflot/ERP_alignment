% Author: David Thinnes
% date:  03/03/2023
% Copyright 2020 by David Thinnes, All rights reserved.

%% set path
close all
clear
clc

method = 'mean';
Fs = 200;
iter_ref = 1;

run('set_path_functions');
load('ERPexample');

[Reg,V] = Var_Alignment_constant(ERPexample,method,iter_ref,Fs);


[synth_ERP, v_ref, reference, groundTruthMap] = generate_synth_ERP();
[Reg_nonconstant,V_nonconstant,initial_denoising, initial_al] = Var_Alignment(synth_ERP,method,iter_ref,Fs,'aniso',false,true);


