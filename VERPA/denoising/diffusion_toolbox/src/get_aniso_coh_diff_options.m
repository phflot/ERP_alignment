% Author   : Philipp Flotho
% Copyright 2019-2020 by Philipp Flotho, All rights reserved.

function options = get_aniso_coh_diff_options()
    options.C = 0.0000001;
    options.alpha = 0.01;
    options.sigma = 0.2;
    options.rho = 4;
    options.tau = 0.24;
    options.iterations = 200;
end

