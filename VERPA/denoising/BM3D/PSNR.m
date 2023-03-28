function PSNRdb = PSNR(x, y, maxval)
    if ~exist('maxval', 'var'), maxval = 255; end
    
    xx=1:size(x,1);
    yy=1:size(x,2);
            
    PSNRdb = zeros(1,size(x,3));
    for ch=1:size(x,3) 
        err = x(xx,yy,ch) - y(xx,yy,ch);
        PSNRdb(ch) = 10 * log10((maxval^2)/mean2(err.^2));    
    end
end