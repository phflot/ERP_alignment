
% Filename : IRANDWT.m
% Date     : 10.01.2018
% Authors  : Ilker Bayram | Manuel C. Kohl

function sig = IRANDWT(cd, ca, N, p, q, r, s, beta)

    % cd : detail coefficients
    % ca : approximation coefficients
    % N : length of the output
    % p, q, r, s : sampling parameters
    % beta : the beta parameter, determines the Q-factor as Q = (2-beta)/beta
    % sig : output signal(s)

    N = N + mod(N,2);
    J = size(cd,2);
    M = size(cd{1}, 1);
    F = RANDWTfilters(N, p, q, r, s, beta, J);
    
    for m = 1:M
        
        X = zeros(1,N);
        w = cell(J+1, 2);
    
        for k = 1:J
            w{k, 1} = real(cd{k}(m, :))';
            w{k, 2} = imag(cd{k}(m, :))';
        end

        w{J+1, 1} = ca(m, :);

        % sampling factors

        PQ = zeros(J,2);
        RS = zeros(J,2);

        for k = 1:J
            p0 = ceil(N*((p/q)^k)); 
            p0 = p0 + mod(p0,2);
            PQ(k,1) = p0;
            RS(k,1) = round(N*((p/q)^(k-1))*r/(2*s));    
        end

        PQ(1,2) = N;
        PQ(2:end,2) = PQ(1:end-1,1);
        RS(1:end,2) = PQ(1:end,2)/2;
        PQ = [1 1;PQ];
        PQk = 1;

        for k = 1:J

            PQk = PQk*PQ(k,1)/PQ(k,2);
            rr = RS(k,1);
            ss = RS(k,2);
            G = conj(F{k,1});
            f = F{k,2};  

            % positive frequencies

            N1 = N*PQk*rr/(2*ss);
            N1 = round(N1);
            d = mod(f{1},N1);
            sub = (fft(w{k,1}-1i*w{k,2},N1))/sqrt(2*N1);
            sub = circshift(sub,-d);
            sub = sub.';    
            dd = length(G);
            X(1+(f{1}:f{1}+dd-1)) = X(1+(f{1}:f{1}+dd-1)) + G.*sub(1:dd);

            % negative frequencies

            g1 = N - f{1};    
            g4 = N - f{1} - dd + 1;     
            sub = (fft(w{k,1} + 1i*w{k,2},N1))/sqrt(2*N1);
            sub = circshift(sub,d-1);
            sub = sub.';    
            X(1+(g4:g1)) = X(1+(g4:g1)) + conj(G(end:-1:1)).*sub(end-dd+1:end);

        end

        sub = fft(w{J+1,1})/sqrt(length(w{J+1,1}));
        H = F{J+1,1};  
        X(1:f{2}) = X(1:f{2}) + sub(1:f{2}).*H(1:f{2});
        X(end-f{2}+2:end) = X(end-f{2}+2:end) + sub(end-f{2}+2:end).*H(f{2}:-1:2);
        x = ifft(X)*sqrt(N);
        sig(m, :) = real(x);
        
    end
    
end
