function ENL = NLEstimate(I,ps,maxiter)
% This function estimate noise level in an Image as specified in paper
% entitled as Single-Image Noise Level Estimation for Blind Denoising by 
% Xinhao Liu, Masayuki Tanaka, and Masatoshi Okutomi.
%
% INPUTS:
%           I = Image
%           ps = Patch Size (Optional); default size in 7
%           maxiter = Number of iteration (Optional); default value is 5
% OUTPUTS:
%           ENL = Estimated Noise Level Can be a single value if I is 
%                 grayscale image or vector of 1 X 3 dimension if I is 
%                 RGB image representing  in red, green and blue channel
%                 
%
% USAGE:
% Estimate Noise Level with default value.
% ENL = NLEstimate(imread('football.jpg'));  % Return Estimated noise level
% for all channel in an image, i.e. for red, green, and blue
% ENL = NLEstimate(rgb2gray(imread('football.jpg')));  % Return Estimated 
% noise level.
%
%
% REFERENCES:
% [1] [Xiang Zhu, and Peyman Milanfar] Automatic Parameter Selection for 
%     Denoising Algorithms Using a No-Reference Measure of Image Content
% [2] [Xinhao Liu Masayuki Tanaka Masatoshi Okutomi] Noise Level Estimation
%     Using Weak Texture Patches of Single Noisy Image 
% [3] [Xinhao Liu, Masayuki Tanaka, and Masatoshi Okutomi] Single-Image 
%     Noise Level Estimation for Blind Denoising
%
% See Also
% http://www.mathworks.com/matlabcentral/fileexchange/36921-noise-level-estimation-from-a-single-image
% 
% Implemented by ASHISH MESHRAM (Meet)
% meetashish85@gmail.com
% Checking input arguments
if nargin<1||isempty(I),error('Missing Input Argument: Specify Image');
else I = im2double(I);end
if nargin<2||isempty(ps),ps = 7;end
if nargin<3||isempty(maxiter),maxiter = 5;end
% Implementation starts here
delta = 0.999; % Confidence Interval, please refer to section III.A from [3]
% Derivative Filter
kh = [-1/2,0,1/2];
% Horizontal Derivative
Ih = imfilter(I,kh,'replicate');
% Cropping to handle borders of each patch correctly
Ih = Ih(:,2:size(Ih,2)-1,:);
% Vertical Derivative
Iv = imfilter(I,kh','replicate');
% Cropping to handle borders of each patch correctly
Iv = Iv(2:size(Iv,1)-1,:,:);
% ------------------- Directional Derivative Operator ------------------- %
% Please refer to reference [3], [1]
Dh = directder(kh,ps,ps);   % Horizontal Direction
Dv = directder(kh',ps,ps);  % Vertical Direction
% ----------------------- Threshold Calculation ------------------------- %
% Evaluating trace of Gradient Matrix
D = trace(Dh'*Dh+Dv'*Dv);
% Intermediate variables 
alfa = (ps*ps)/2;
beta = D/alfa;
% Computing threshold from Inverse gamma function
tau0 = gaminv(delta,alfa,beta);
% Preallcoating Noise Level array for three iamge channel
ENL = zeros(1,size(I,3));
sigma = zeros(1,maxiter);
for chan = 1:size(I,3)
    % Selecting particular channel (from red, green, and blue) of image,
    % and converting entire image into vector of specified patch size
    X = im2col(I(:,:,chan),[ps ps]);
    Xh = im2col(Ih(:,:,chan),[ps ps-2]);
	Xv = im2col(Iv(:,:,chan),[ps-2 ps]);
    Xtr = sum(vertcat(Xh,Xv));
    
    %---------------- Initial Noise Level Estimation ---------------------%
    % Computing Covariance Matrix
    C = (X*X')/(size(X,2)-1);
    % Compute Eigen Value of covariance matrix
    EV = eig(C);
    % Initial variance of Image
    sigma(1) = EV(1);
    
    %----------------- Iterative Noise Level Estimation ------------------%
    for k = 2:maxiter
        % Updating threshold equation 18 of [3]
        tau = sigma(k-1)*tau0;
        % Get Weak Texture Patch Index Using Selection Criteria, please
        % refer to section III.A equation no. 17 from [3]
        p = (Xtr<tau);
        Xtr = Xtr(:,p);
        % Selecting Weak Texture Patch
        WTP = X(:,p);
        
        %----- Noise Level Estimation for Selected Weak Texture Patch ----%
        % Covariance Matrix
        CWTP = (WTP*WTP')/(size(WTP,2)-1);
        % Compute Eigen Value of Covariance Matrix
        EVWTP = eig(CWTP);
        % Selecting minimum eigen value as variance of image and updating
        % initial value
        sigma(k) = EVWTP(1);
        
        % Specifying condtion for stoping loop
        if abs(sigma(k) - sigma(k-1)) <= 0.0001
            sig = sigma(k);
            break;
        end
    end
    % Estimated Noise Level
    ENL(chan) = sqrt(sig);
end