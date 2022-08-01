% Author   : Philipp Flotho
% Copyright 2019-2020 by Philipp Flotho, All rights reserved.

function filtered = aniso_diff_filt(img, options)
    global use_gpu;       
    try 
        gpuArray(1);
        use_gpu = true;
        fprintf('Using gpu\n');
    catch
        use_gpu = false;
        fprintf('Not using gpu\n');
    end

    filtered = zeros(size(img), 'double');
    
    tic
    for i = 1:size(img, 3)
        filtered(:, :, i) = ...
            diffusion_wrapper(im2double(img(:, :, i)), options);
    end
    toc
end

function filtered = diffusion_wrapper(img, options)
    global use_gpu;

    [m, n] = size(img);
    
    if (use_gpu)
        u = gpuArray(zeros(m, n, 2));

        u(:, :, 1) = gpuArray(img);
    else
        u = zeros(m, n, 2);

        u(:, :, 1) = img;
    end
    
    switch options.diffusivity
        case '1'
            diffusivity = @diffusivity1;
        case '2'
            diffusivity = @diffusivity2;
        case 'weickert'
            diffusivity = @diffusivity3;
    end

    for i = 1:options.iterations
        [a, b, c] = get_diff_tensor(u(:, :, 1), ...
            options.sigma, options.K, diffusivity);
        u = diffusion_step(u, a, b, c, options.tau);
    end
    filtered = gather(u(:, :, end));
end

function u = diffusion_step(u, a, b, c, tau)
    
    global use_gpu;

    t = 2;
    
    [m, n, ~] = size(u);
    
    if (use_gpu)
        c1 = gpuArray(2:m-1);
        c2 = gpuArray(2:n-1);
        l = gpuArray(1:n-2);
        to = gpuArray(1:m-2);
        r = gpuArray(3:n);
        bo = gpuArray(3:m);
    else
        c1 = 2:m-1;
        c2 = 2:n-1;
        l = 1:n-2;
        to = 1:m-2;
        r = 3:n;
        bo = 3:m;
    end

    u_tl = u(to, l, t-1) .* 1/4 .* (...
          abs(b(to, l)) + b(to, l)...
        + abs(b(c1, c2)) + b(c1, c2));

    u_tr = u(to, r, t-1) .* 1/4 .* (...
          abs(b(to, r)) - b(to, r)...
        + abs(b(c1, c2)) - b(c1, c2));

    u_bl = u(bo, l, t-1) .* 1/4 .* (...
          abs(b(bo, l)) - b(bo, l)...
        + abs(b(c1, c2)) - b(c1, c2));

    u_br = u(bo, r, t-1) .* 1/4 .* (...
          abs(b(bo, r)) + b(bo, r)...
        + abs(b(c1, c2)) + b(c1, c2));

    u_tc = u(to, c2, t-1) .* 1/2 .* (...
          c(to, c2) + c(c1, c2)...
        -(abs(b(to, c2)) + abs(b(c1, c2))));

    u_bc = u(bo, c2, t-1) .* 1/2 .* (...
          c(bo, c2) + c(c1, c2)...
        -(abs(b(bo, c2)) + abs(b(c1, c2))));

    u_cl = u(c1, l, t-1) .* 1/2 .* (...
          a(c1, l) + a(c1, c2)...
        -(abs(b(c1, l)) + abs(b(c1, c2))));

    u_cr = u(c1, r, t-1) .* 1/2 .* (...
          a(c1, r) + a(c1, c2)...
        -(abs(b(c1, r)) + abs(b(c1, c2))));

    u_cc = u(c1, c2, t-1) .* (1/2 .* (...
        -(a(c1, l) + 2*a(c1, c2) + a(c1, r))...
        + abs(b(c1, l)) + abs(b(c1, r)) + ...
            abs(b(to, c2)) + abs(b(bo, c2)) + ...
            2 .* abs(b(c1, c2))...
        - (c(to, c2) + 2.*c(c1, c2) + c(bo, c2))) ...
        + 1/4.*(...
        - (abs(b(bo, l)) - b(bo, l) + ...
            abs(b(bo, r)) + b(bo, r))...
        - (abs(b(to, l)) + b(to, l) + ...
            abs(b(to, r)) - b(to, r))));

    u(2:end-1, 2:end-1, t) = (...
          u_tl + u_tc + u_tr + u_cl + u_cc + u_cr + u_bl + u_bc + u_br ...
          ) .* tau + u(2:end-1, 2:end-1, t-1);
    u(:, :, t-1) = u(:, :, t);
end

function [a, b, c] = get_diff_tensor(img, sigma, K, diffusivity)

    [dx, dy] = gradient(imgaussfilt(img, sigma));
    dy = 0.15 * dy;
    grad_mag = sqrt(dx.^2 + dy.^2);

    v1_1 = dx ./ (grad_mag + 0.000001);
    v1_2 = dy ./ (grad_mag + 0.000001);

    v2_1 = dy ./ (grad_mag + 0.000001);
    v2_2 = -dx ./ (grad_mag + 0.000001);

    lambda1 = diffusivity(img, K);
    lambda2 = 1;

    a = lambda1 .* v1_1 .* v1_1 + ...
        lambda2 .* v2_1 .* v2_1;
    b = lambda1 .* v1_1 .* v1_2 + ...
        lambda2 .* v2_1 .* v2_2;
    c = lambda1 .* v1_2 .* v1_2 + ...
        lambda2 .* v2_2 .* v2_2;
end

function d = diffusivity1(img, K)
% K is the edge parameter
    [grad_mag, ~] = imgradient(img, 'central');
    d = exp(-grad_mag./K);
end

function d = diffusivity2(img, K)
% K is the edge parameter
    [grad_mag, ~] = imgradient(img, 'central');
    d = 1 ./ (1 + (grad_mag./K).^2);
end

function d = diffusivity3(img, K)
% weickert diffusivity
    global use_gpu;

    [grad_mag, ~] = imgradient(img, 'central');

    if use_gpu
        d = gpuArray(ones(size(img)));
    else
        d = ones(size(img));
    end
    
    idx = grad_mag.^2 ~= 0;
    d(idx) = 1 - exp(...
            (-3.31488/255) ./ (grad_mag(idx)./K).^8);
end