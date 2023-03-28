% Filename : unidirectionalGaussianMeans.m
% Date     : 18.01.2018
% Author   : Manuel C. Kohl

function filteredERPimage = unidirectionalGaussianMeans(ERPimage, sigma)

    % Padding ERP image with mirrored copies

    [nResponses, nSamples] = size(ERPimage);
    ERPimage = mirrorPadding(ERPimage, [sigma 0]);
    nPaddedResponses = size(ERPimage, 1);

    % Convolving with Gaussian window across responses

    sweeps = 1:nPaddedResponses;
    window = exp(-((sweeps-mean(sweeps))/sigma).^2);
    window = window / sum(window);
    for k = 1:nSamples
        ERPimage(:, k) = conv(ERPimage(:, k), window, 'same');
    end
	
	% Removing padding

    filteredERPimage = ERPimage(sigma+1:sigma+nResponses, :);

end
