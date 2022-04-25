
% Filename : ASPRdenoising.m
% Date     : 02.02.2018
% Authors  : Manuel C. Kohl

function denoisedERPimage = ASPRdenoising(ERPimage)

	nSamples = size(ERPimage, 2);

	% Parameters

	p = 7;
	q = 9;
	D = p/q;
	r = 1;
	s = 2;
    % R = r/(s*(1-D));
	Q = 3;
	beta = 2/(Q+1);
	J = floor(log(4*s/(r*(nSamples+mod(nSamples, 2))))/log(D));
	sigma = 5;
	mu = 1;
	it = 10;
	lambda = 0.01;
	
	% Performing RANDWT
            
	[cd, ca] = RANDWT(ERPimage, p, q, r, s, beta, J);
	cdf = cd;
	caf = ca;
	
	% Smoothing approximation coefficients
	
	caf = unidirectionalGaussianMeans(caf, sigma);
	
	% Calculating localized phase stability of detail coefficients
	
	lps = cell(1, J);
	
	for k = 1:J
		lps{k} = abs(unidirectionalGaussianMeans(exp(1i*angle(cd{k})), sigma));
	end
	
	% Performing amplitude and phase regularization
	
	for k = 1:J
		cam = abs(cdf{k});
		cph = exp(1i*angle(cdf{k}));
		cam = cam .* lps{k};
		cam = unidirectionalGaussianMeans(cam, sigma);
		cph = isotropicPhaseRegularization(cph, mu, it, lambda);
		cdf{k} = cam .* exp(1i*angle(cph));
	end
	
	% Performing inverse RANDWT
	
	denoisedERPimage = IRANDWT(cdf, caf, nSamples, p, q, r, s, beta);
    denoisedERPimage = denoisedERPimage(:, 1:nSamples);
	
end
