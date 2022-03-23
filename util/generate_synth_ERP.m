% Date     : 17.07.2020
% Author   : Philipp Flotho

function [img, v_ref, reference, groundTruthMap] = generate_synth_ERP(rng_idx, v_ref)
    run('../set_path.m');

    reference = [0 0 0 0 0 0 0 0 0 0 0 0 0 0.0073 0.0286 0.0629 0.108 0.1613 0.2199 0.2801 0.3387 0.392 0.4226 0.4141 0.367 0.2864 0.1796 0.0532 -0.0851 -0.2262 -0.3611 -0.4806 -0.5764 -0.6413 -0.6698 -0.6583 -0.6054 -0.5121 -0.3816 -0.2263 -0.0603 0.1093 0.2748 0.4292 0.5656 0.6784 0.7633 0.8172 0.8536 0.8865 0.9157 0.941 0.9619 0.9785 0.9904 0.9976 1 0.9976 0.9904 0.9785 0.9619 0.941 0.9157 0.8865 0.8536 0.8172 0.7778 0.7357 0.6913 0.6451 0.5975 0.549 0.5 0.451 0.4025 0.3549 0.3087 0.2643 0.2222 0.1828 0.1464 0.1135 0.0843 0.059 0.0381 0.0215 0.0096 0.0024 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

    [~, idx] = findpeaks(reference);
    
    kernel_size = 71;

    n_trials = 400;
    width = length(reference);

    tmp = sin(linspace(0, 1, n_trials) * 3 * pi)-1;

    if nargin >= 1
        rng(rng_idx);
    else
        rng(20);
    end
    
    img = repmat(reference, n_trials, 1);
    
    if nargin < 2
        v_ref = zeros(size(img));
        
        % global displacements: 
        max_disp = 1;
        v_ref = zeros(size(img));
        prev_rand = 0;
        for i = 1:n_trials
            prev_rand = (rand - 0.5) * max_disp + prev_rand;
            v_ref(i, :) = prev_rand;
        end

        v_ref(1, :) = 0;
        
        % local displacements
        for i = 1:n_trials
            v_ref(i, idx(1)-10:idx(1)+ 10) = v_ref(i, 1) + normrnd(0, 0.5);
            v_ref(i, idx(2)-15:idx(2)+ 15) = v_ref(i, 1) + normrnd(0, 2);
        end
        v_ref = imgaussfilt2(v_ref, [0, 5]);
    end
    
    m = width;
    n = n_trials;
    [X, Y] = meshgrid(1:m, 1:n);
    F = scatteredInterpolant(...
        X(:) + v_ref(:), ...
        Y(:) + zeros(m * n, 1), img(:));
    groundTruthMap = F(X, Y);
    img = groundTruthMap + ...
        imgaussfilt2(normrnd(0, 2, n_trials, width), [0, 2]) + ...
        normrnd(0, 0.5, n_trials, width);
end

