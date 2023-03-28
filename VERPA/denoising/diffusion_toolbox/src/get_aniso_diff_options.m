% Author   : Philipp Flotho
% Copyright 2019-2020 by Philipp Flotho, All rights reserved.

function options = get_aniso_diff_options()
    options.sigma = 2;
    options.K =  0.01;
    options.tau = 0.24;
    options.iterations = 200;
    options.diffusivity = '1';
end