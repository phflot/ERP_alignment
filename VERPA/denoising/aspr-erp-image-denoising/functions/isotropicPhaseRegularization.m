
% Filename : isotropicPhaseRegularization.m
% Date     : 18.01.2018
% Authors  : Daniel J. Strauss | Manuel C. Kohl

function dataRegIso = isotropicPhaseRegularization(data, mu, it, lambda)

    [n, m] = size(data);
    dataRegIso = data;
    abl_iso = zeros(size(data));
    data_mirror_iso = mirrorPadding(dataRegIso, [1 1]);

    for w = 1:it
        f0_mirror_iso = mirrorPadding(dataRegIso, [1 1]);
        A = 2+16*mu;
        for jj = 2:n+1
            for ii = 2:m+1
                b = 2*data_mirror_iso(jj,ii)+2*mu*(f0_mirror_iso(jj,ii-1)+f0_mirror_iso(jj,ii+1)...
                    +f0_mirror_iso(jj-1,ii)+f0_mirror_iso(jj+1,ii)+f0_mirror_iso(jj-1,ii-1)...
                    +f0_mirror_iso(jj+1,ii+1)+f0_mirror_iso(jj-1,ii+1)+f0_mirror_iso(jj+1,ii-1));
                dataRegIso(jj-1,ii-1) = A\b;
                abl_iso(jj-1,ii-1) = 2*(f0_mirror_iso(jj,ii)-data_mirror_iso(jj,ii))...
                    -2*mu*(f0_mirror_iso(jj,ii-1)+f0_mirror_iso(jj,ii+1)...
                    -2*f0_mirror_iso(jj,ii))...
                    -2*mu*(f0_mirror_iso(jj-1,ii)+f0_mirror_iso(jj+1,ii)...
                    -2*f0_mirror_iso(jj,ii))...
                    -2*mu*(f0_mirror_iso(jj-1,ii-1)+f0_mirror_iso(jj+1,ii+1)...
                    -2*f0_mirror_iso(jj,ii))...
                    -2*mu*(f0_mirror_iso(jj-1,ii+1)+f0_mirror_iso(jj+1,ii-1)...
                    -2*f0_mirror_iso(jj,ii));
            end
        end
        dataRegIso = dataRegIso - lambda.*abl_iso;
    end

end
