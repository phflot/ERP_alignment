function D = directder(H, rps, cps)
% This function creates Directional derivative as specified in paper
% entitlled as Single-Image Noise Level Estimation for Blind Denoising by 
% Xinhao Liu, Masayuki Tanaka, and Masatoshi Okutomi.
%
% INPUT:    H = 
%           rps = 
%           cps = 
% OUTPUT:
%           D = 
%
% USAGE:
% H = [1 2 3];
% 
% D = directder(H);
% D = directder(H,3);
% D = directder(H,3,5);
% REFERENCE:
% [1] Single-Image Noise Level Estimation for Blind Denoising Xinhao Liu, 
%     Masayuki Tanaka, and Masatoshi Okutomi
%
% Implemented by ASHISH MESHRAM (Meet)
% meetashish85@gmail.com
% Checking input arguments
if nargin<1||isempty(H),error('Missing Input Argument');end
if nargin<2||isempty(rps),rps = 7;cps = 7;end
if nargin<3||isempty(cps),cps = rps;end
% Implementation starts here {Please refer to reference)
V = zeros(1,rps*cps);
V(1:length(H)) = H;
D = toeplitz(V);
