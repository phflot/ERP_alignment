% Author   : Philipp Flotho
% Copyright 2019-2020 by Philipp Flotho, All rights reserved.

function filtered = perona_malik_diff(img, options)
% perona Malik diffusion
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
    
    tau = options.tau;
    iterations = options.iterations;
    hs = 1;

    [m, n] = size(img);
    if (use_gpu)
        u = gpuArray(zeros(m, n, iterations));

        u(:, :, 1) = gpuArray(img);
    else
        u = zeros(m, n, iterations);
        u(:, :, 1) = img;
    end
       
    u(:, :, 1) = gpuArray(img);

    for t = 2:iterations
        switch options.diffusivity
            case 'weickert'
                d = diffusivity3(imgaussfilt(...
                    u(:, :, t-1), options.sigma), options.K);
            case '1'
                d = diffusivity1(imgaussfilt(...
                    u(:, :, t-1), options.sigma), options.K);
            otherwise
                d = diffusivity2(imgaussfilt(...
                    u(:, :, t-1), options.sigma), options.K);
        end
        [dx, dy] = imgradientxy(d, 'central');
        
        [ux, uy] = imgradientxy(u(:, :, t-1), 'central');

        uxx = hs * (u(3:end, 2:(end - 1), t - 1) ...
                  - 2 * u(2:(end - 1), 2:(end - 1), t - 1) ...
                  + u(1:(end - 2), 2:(end - 1), t - 1));
        uyy = hs * (u(2:(end - 1), 3:end, t - 1) ...
                  - 2 * u(2:(end - 1), 2:(end - 1), t - 1) ...
                  + u(2:(end - 1), 1:(end - 2), t - 1));    

        u(2:end-1, 2:end-1, t) = tau .* (...
                d(2:end-1, 2:end-1) .* uxx + dx(2:end-1, 2:end-1) .* ux(2:end-1, 2:end-1)...
              + d(2:end-1, 2:end-1) .* uyy + dy(2:end-1, 2:end-1) .* uy(2:end-1, 2:end-1))...
              + u(2:end-1, 2:end-1, t-1);
    end

    u(isnan(u(:))) = 0;
    if use_gpu
        filtered = gather(mat2gray(u(:, :, end)));
    else
        filtered = u(:, :, end);
    end
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
    [grad_mag, ~] = imgradient(img, 'central');

    if (grad_mag^2 == 0)
        d = 1;
    else
        d = 1 - exp(...
            (-3.31488/255) ./ (grad_mag/K).^8);
    end
end
