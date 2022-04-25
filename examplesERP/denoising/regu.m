function f01_iso = regu(data,mu,it,lambda)

[n m] = size(data);

f01_iso = data;
abl_iso = zeros(size(data));

f01_vert = data;
abl_vert = zeros(size(data));

data_mirror_iso = mirror_xy(f01_iso,1,1); 
data_mirror_vert = mirror_xy(f01_vert,1,1); 

for w = 1:it % iteration
    f0_mirror_iso = mirror_xy(f01_iso,1,1);% f0
    f0_mirror_vert = mirror_xy(f01_vert,1,1);
    A = 2+16*mu;
    for jj = 2:n+1
        for ii = 2:m+1
                b = 2*data_mirror_iso(jj,ii)+2*mu*(f0_mirror_iso(jj,ii-1)+f0_mirror_iso(jj,ii+1)...
                +f0_mirror_iso(jj-1,ii)+f0_mirror_iso(jj+1,ii)+f0_mirror_iso(jj-1,ii-1)...
                +f0_mirror_iso(jj+1,ii+1)+f0_mirror_iso(jj-1,ii+1)+f0_mirror_iso(jj+1,ii-1));
                 f01_iso(jj-1,ii-1) = A\b; 
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
    f01_iso = f01_iso-lambda.*abl_iso;
end



% % plot
%  figure;
%  subplot(1,4,1)
%  imagesc(data);
%  title('noisy');
%  ylabel('Number of Sweeps');
%  xlabel('Samples');
%  colorbar;
%  
% 
%  
%  subplot(1,4,2)
%  plot(data(25,:));
%  hold on;
%  plot(data(100,:));
%  title('sweeps 25 and 100');
%  xlabel('Samples');
% 
% 
% 
%  
%  subplot(1,4,3)
%  imagesc(f01_iso);
%  title('Nonlinear Diffusion');
%  xlabel('Samples');
%  colorbar;
% 
%  
%  subplot(1,4,4)
%  plot(f01_iso(25,:));
%  hold on;
%  plot(f01_iso(100,:));
%  title('sweeps 25 and 100');
%  xlabel('Samples');
%  colormap(hot); 





function BB = mirror_xy(B,Llo,Lru)
%------------------------------------------------------------------------
% creates boundaries by mirroring
% Llo: number of pixels added at the top and on the left side
% Lru: number of pixels added at the bottom and on the right side

[M,N] = size(B);

%Llo,Lru <= M,N
%if max(Llo,Lru) > min(M,N)
%    error('Filter length exceeds image size.') 
%end    

%index vectors
mind = [fliplr(1:Llo) (1:M) fliplr(M-Lru+1:M)];
nind = [fliplr(1:Llo) (1:N) fliplr(N-Lru+1:N)];

BB = B(mind,nind);
% BB = B(mind,1:N);
