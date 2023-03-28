function [PMSE, z_final] = kernel_CG_W (y,y_est, ...
    Hf, W,  gama, eta, numiter,  zh)

%% MATLAB function for SD implementation of kernel based denoising algorithm
[N,M] = size(y);
z_est = y_est;
bnn1 = W*y(:); bnn = reshape(bnn1,N,M);
b = real(ifft2(fft2((1+gama)*y-gama*bnn).*conj(Hf)));
bb1 = W*z_est(:); b11 = reshape(bb1,N,M);
z1 = real(ifft2(fft2(z_est).*Hf));
bnn1 = W*z1(:); bnn = reshape(bnn1,N,M);
cc = real(ifft2(fft2(-gama*bnn+(1+gama)*z1).*conj(Hf)))+eta*(z_est-b11);
r = b-cc; rho = norm(r,'fro')^2;
p = zeros(N,M); w1 = zeros(N,M); zh1 = zeros(N,M); z_est1 = zeros(N,M);

for i = 1:numiter
    if i==1
        p=r; rho1=rho; p=r+(rho/rho1)*p;
    z1 = real(ifft2(fft2(p).*Hf));

        bb1 = W*p(:);b11 = reshape(bb1,N,M);
        bnn1 = W*z1(:); bnn = reshape(bnn1,N,M);
w1 = real(ifft2(fft2(-gama*bnn+(1+gama)*z1).*conj(Hf)))+eta*(p-b11);
        alpha1 = rho/sum(sum(p.*w1));
        z_est = z_est+alpha1*p;
        r = r-alpha1*w1;
        rho1 = rho;
%         rho = norm(r(:))^2;
rho = norm(r,'fro')^2;
        zh1 = zh - real(ifft2(fft2(z_est).*Hf));
        PMSE(i) = norm(zh1, 'fro');
    else
        p=r+(rho/rho1)*p;
    z1 = real(ifft2(fft2(p).*Hf));
        bb1 = W*p(:);
        b11 = reshape(bb1,N,M);
        bnn1 = W*z1(:);
        bnn = reshape(bnn1,N,M);

    w1 = real(ifft2(fft2(-gama*bnn+(1+gama)*z1).*conj(Hf)))+eta*(p-b11);

        alpha1 = rho/sum(sum(p.*w1));
        z_est1 = z_est;
        z_est = z_est+alpha1*p;
        r = r-alpha1*w1;
        rho1 = rho;
        rho = norm(r(:))^2;
       zh1 = zh - real(ifft2(fft2(z_est).*Hf));
        PMSE(i) = norm(zh1, 'fro');
        if(PMSE(i)>PMSE(i-1))
            break;
        end
    end
end
z_final = z_est1;