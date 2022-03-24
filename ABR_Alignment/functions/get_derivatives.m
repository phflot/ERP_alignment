% Date     : 17.07.2019
% Author   : Philipp Flotho

function [fy, fyy, f_ref_y] = get_derivatives(f, f_ref, hx, hy)
    
    f = f - min(f, [], 1);
    f = f./ max(f, [], 1);
    
    f_ref = f_ref - min(f_ref(:));
    f_ref = f_ref./ max(f_ref(:));

    f = set_boundary(padarray(f, [1 1], 'symmetric'));

    fyy = zeros(size(f));
    
    fyy(2:end-1, 2:end-1) = ...
        (f(1:end-2, 2:end-1) - 2 * f(2:end-1, 2:end-1) ...
        + f(3:end, 2:end-1)) ./ hy^2;
    
    [~, fy] = gradient(set_boundary(f), hx, hy);
    f_ref_y = gradient([f_ref(2); f_ref; f_ref(end-1)], hy);
end

function f = set_boundary(f)
    f(:, 1) = f(:, 3);
    f(:, end) = f(:, end - 2);
    f(1, :) = f(3, :);
    f(end, :) = f(end - 2, :);
end