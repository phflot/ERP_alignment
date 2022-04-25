
% Filename : RANDWTfilters.m
% Date     : 10.01.2018
% Authors  : Ilker Bayram | Manuel C. Kohl

function F = RANDWTfilters(N, p, q, r, s, beta, J)

    % N : input signal length
    % p, q, r, s : sampling parameters
    % beta : the beta parameter, determines the Q-factor as Q = (2-beta)/beta
    % J : number of subbands

    N = N + mod(N,2);

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
    Hk = ones(1,N);

    for k = 1:J
        
        PQk = PQk*PQ(k,1)/PQ(k,2);
        pp = PQ(k+1,1);
        qq = PQ(k+1,2);
        rr = RS(k,1);
        epsi = (N/64)*(pp-qq + beta*qq)/(pp+qq);
        f{1} = ((1-beta)*N/2 + epsi);
        f{2} = (N/2*pp/qq);
        f{3} = (N/2 - epsi);
        f{4} = (N/2 + epsi);

        for n = 1:4
            f{n} = (f{n}*PQk);
        end
        
        [H, G] = MakeFilters(f);
        f{1} = ceil(f{1});  
        f{2} = floor(f{2});
        f{4} = floor(f{4});

        % filter for positive frequencies
        
        Gk = Hk(1+(f{1}:f{4})).*G(1:end);
        dd = min(rr,length(Gk));
        F{k,1} = Gk(1:dd);
        F{k,2} = f;

        % update lowpass filter
        
        Hk(1+(f{1}:f{2})) = Hk(1+(f{1}:f{2})).*H(1+(f{1}:f{2})); 
        Hk(1+(f{2}+1:f{4})) = 0;

    end
    
    F{J+1,1} = Hk(1:f{2});
    
end

function [H, G] = MakeFilters(f)

    % MAKE H0 and G0
    
    w = (0:floor(f{2}));
    k_pass = (w < f{1});      
    k_trans = (w >= f{1});
    b = (f{2}-f{1})/pi;
    w_scaled = (w - f{1})/b;
    H = zeros(size(w));
    H(k_pass) = 1;
    H(k_trans) = (1+cos(w_scaled(k_trans))) .* sqrt(2-cos(w_scaled(k_trans)))/2;
    seq = sqrt(1 - H(k_trans).^2);
    w = (ceil(f{1}):floor(f{4}));
    k_pass = (w <= f{3}) & (w >= f{2});
    k_trans1 = (w <= f{2});
    k_trans2 = (w >= f{3});
    b = (f{4}-f{3})/pi;
    
    if b > 0
        w_scaled = (w - f{3})/b;
    else
        w_scaled = 0*w;
    end

    G = zeros(size(w));
    G(k_pass) = 1;
    G(k_trans1) = seq;
    G(k_trans2) = (1+cos(w_scaled(k_trans2))) .* sqrt(2-cos(w_scaled(k_trans2)))/2;
    
end
