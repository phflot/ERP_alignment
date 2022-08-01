% Author   : Philipp Flotho
% Copyright 2019-2021 by Philipp Flotho, All rights reserved.

function filtered = aniso_coh_diff_filt(img, options)
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
    for i = 1:options.iterations
        [a, b, c] = get_diff_tensor(u(:, :, 1), ...
            options.sigma, options.C, options.alpha, options.rho);
        u = diffusion_step(u, a, b, c, options.tau);
    end
    if use_gpu
        filtered = gather(u(:, :, end));
    else
        filtered = u(:, :, end)
    end
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

function [a, b, c] = get_diff_tensor(img, sigma, C, alpha, rho)
    global use_gpu;
    [dx, dy] = gradient(imgaussfilt(img, sigma));
%     dy = dy * 0.15;

    [m, n, ~] = size(img);
    
    % structure tensor: 
    a = imgaussfilt(dx .* dx, rho);
    b = imgaussfilt(dx .* dy, rho);
    c = imgaussfilt(dy .* dy, rho);

    mu1 = (a + c + ...
        sqrt((c-a).^2 + 4 .* b .*b))./2;
    mu2 = (a + c - ...
        sqrt((c-a).^2 + 4 .* b .*b))./2;

    v1_1 = 2 .* b;
    v1_2 = c - a + ...
        sqrt((c-a).^2 + 4 .* b .*b);
    mag = sqrt(v1_1.^2 + v1_2.^2);

    v1_1 = v1_1 ./ mag;
    v1_2 = v1_2 ./ mag;
    v1_1(mag < 0.000001) = 1;
    v1_2(mag < 0.000001) = 0;

    if (use_gpu)
        lambda1 = alpha .* gpuArray(ones(m, n));
        lambda2 = alpha .* gpuArray(ones(m, n));
    else
        lambda1 = alpha .* ones(m, n);
        lambda2 = alpha .* ones(m, n);
    end
        
    idx = mu1 ~= mu2;
    lambda2(idx) = ...
        alpha + (1 - alpha) .* ...
        exp(-C./((mu1(idx) - mu2(idx)).*(mu1(idx) - mu2(idx))));

    a = lambda1 .* v1_1 .* v1_1 + ...
        lambda2 .* v1_2 .* v1_2;
    c = lambda1 + lambda2 - a;
    b = (lambda1 - lambda2) .* v1_1 .* v1_2;
end