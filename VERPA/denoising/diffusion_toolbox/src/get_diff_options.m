% Author   : Philipp Flotho
% Copyright 2019-2020 by Philipp Flotho, All rights reserved.

function options = get_diff_options()
    options.K = 0.01; 
    options.sigma = 2;
    options.iterations = 50;
    options.tau = 0.05;
    options.diffusivity = '1';
end

