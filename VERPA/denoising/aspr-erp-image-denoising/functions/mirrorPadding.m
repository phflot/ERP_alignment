
% Filename : mirrorPadding.m
% Date     : 18.01.2018
% Author   : Daniel J. Strauss | Manuel C. Kohl

function paddedERPimage = mirrorPadding(ERPimage, paddingSize)

    [m, n] = size(ERPimage);
    mIndices = [fliplr(1:paddingSize(1)) (1:m) fliplr(m-paddingSize(1)+1:m)];
    nIndices = [fliplr(1:paddingSize(2)) (1:n) fliplr(n-paddingSize(2)+1:n)];
    paddedERPimage = ERPimage(mIndices, nIndices);

end
