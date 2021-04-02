function omega = initialize_omega(init,N,K)


omega = zeros(N, 2, K);

% Initialization of omega_k
switch init
    case 0

        for k = 1:K
            omega(1,1,k) = 0.25*cos(pi*(k-1)/K);
            omega(1,2,k) = 0.25*sin(pi*(k-1)/K);
        end
        
        % Case 1: random on half-plane
    case 1
        for k=1:K
            omega(1,1,k) = rand()-1/2;
            omega(1,2,k) = rand()/2;
        end

end