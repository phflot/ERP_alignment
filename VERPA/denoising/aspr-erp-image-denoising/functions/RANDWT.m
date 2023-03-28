
% Filename : RANDWT.m
% Date     : 10.01.2018
% Authors  : Ilker Bayram | Manuel C. Kohl

function [cd, ca] = RANDWT(sig, p, q, r, s, beta, J)

    % sig : input signal(s)
    % p, q, r, s : sampling parameters
    % beta : the beta parameter, determines the Q-factor as Q = (2-beta)/beta
    % J : number of subbands
    % cd : detail coefficients
    % ca : approximation coefficients

    M = size(sig, 1);
    
    for m = 1:M
    
        x = sig(m, :);
        L = length(x); 
        N = L + mod(L,2);
        x = [x zeros(1,N-L)];

        if (N * ((p/q)^J))*r/(2*s) < 2
            error('Too many levels - reduce J');
        end

        F = RANDWTfilters(N, p, q, r, s, beta, J);
        X = fft(x)/sqrt(N);

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
            G = F{k,1};
            f = F{k,2};

            % positive frequencies

            dd = length(G);
            sub = G(1:end).*X(1+(f{1}:f{1}+dd-1));
            N1 = N*PQk*rr/(2*ss);
            N1 = round(N1);
            sub2 = zeros(1,N1);    
            sub2(1:length(sub)) = sub;
            d = mod(f{1},N1);    
            sub2 = circshift(sub2.',d);
            sub2 = ifft(sub2);

            % negative frequencies

            g1 = N - f{1};    
            g4 = N - f{1} - dd + 1;
            sub = conj(G(end:-1:1)).*X(1+(g4:g1));
            sub3 = zeros(1,N1);
            sub3(end-length(sub)+1:end) = sub;
            sub3 = circshift(sub3.',-(d-1));
            sub3 = ifft(sub3);
            w{k,1} = (sub2 + sub3)*sqrt(N1/2);
            w{k,2} = 1i*(sub2 - sub3)*sqrt(N1/2);

        end

        H = F{J+1};
        N1 = PQ(J+1,1); 
        sub = [X(1+(0:f{2}-1)).*H(1+(0:f{2}-1)) 0 X(end-f{2}+2:end).*H(f{2}:-1:2)];
        sub2 = ifft(sub)*sqrt(N1);
        w{J+1,1} = (sub2);   

        for k = 1:J
            cd{k}(m, :) = w{k, 1} + 1i*w{k, 2};
        end

        ca(m, :) = w{J+1, 1};
        
    end
           
end
