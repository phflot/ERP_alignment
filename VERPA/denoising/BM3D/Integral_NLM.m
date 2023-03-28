function [Kernel] = Integral_NLM(InImage,krad,wrad, h)
% Function Integral_NLM_sparse : Computes the kernel similarity matrix
% using Integram Images
% Inputs:
%  InImage = Input Image
%  krad = The patch radius  for computation of the weights
%  wrad = search neighborhood radius for computaion of the weights
%  h = smoothing (scaling) parmeter for computation of the weights

% Outputs:
%  Kernel = Sparse kernel Similarity Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,M] = size(InImage);
indmat = reshape(1:N*M,N,M);
wsz = 2*wrad+1; ihwl = 2*krad+1; iwl = ihwl*ihwl;
Kernel1 = zeros(N*M*wsz^2,1);
indx = repmat((1:N*M)',wsz^2,1);
% Padding the input image and the index image
PaddedImg = padarray(InImage,[krad+wrad,krad+wrad],'symmetric');
imgind = padarray(indmat, [krad+wrad,krad+wrad]);
indy = zeros(N*M*wsz^2,1);
q = 0;
for xShift = -wrad:wrad
    for yShift = -wrad:wrad

        ssd = IntegralImage(PaddedImg,xShift,yShift); % integral image
        % deriving the patch squared difference between the corresonding
        % pairs of pixels
        PatchDist = ssd(wrad+2*krad+1:end-wrad,wrad+2*krad+1:end-wrad)+...
            ssd(wrad:end-2*krad-wrad-1,wrad:end-2*krad-wrad-1)-...
            ssd(wrad:end-2*krad-wrad-1,wrad+2*krad+1:end-wrad)-...
            ssd(wrad+2*krad+1:end-wrad,wrad:end-2*krad-wrad-1);
        % Computing the similarity weights
        w = exp(-PatchDist/(h^2*iwl));
        % Modifying the weights for the border pixels
        vind = imgind((krad+wrad+1+xShift):(krad+wrad+xShift+N),(krad+wrad+1+yShift):(krad+wrad+yShift+M));
        [ff] = find(vind==0);w(ff)=0;vind(ff)=1;
        indy (q*N*M+1:q*N*M+N*M)= vind(:);
        Kernel1 (q*N*M+1:q*N*M+N*M) = w(:);
        q=q+1;
    end
end

Kernel = sparse(indx, indy, Kernel1, (N*M), (N*M),numel(Kernel1));

end



