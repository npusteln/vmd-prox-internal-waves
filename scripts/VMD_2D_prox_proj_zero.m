function [u, u_hat, omega,crit] = VMD_2D_prox_proj_zero(signal, Alpha, rho,eta, beta, K, omega, tol,N)


% Resolution of image
[Hy,Hx] = size(signal);
[X,Y] = meshgrid((1:Hx)/Hx, (1:Hy)/Hy);


% Spectral Domain discretization
fx = 1/Hx;
fy = 1/Hy;
freqs_1 = X - 0.5 - fx;
freqs_2 = Y - 0.5 - fy;





% Construct f and f_hat
f_hat = fftshift(fft2(signal));

% Storage matrices for (Fourier) modes. All iterations are not recorded.
u_hat = zeros(Hy,Hx,K);
u_hat_old = u_hat;

% Storage matrices for (Fourier) Lagrange multiplier.
mu_hat = zeros(Hy,Hx);




%% Main loop for iterative updates

% Stopping criteria tolerances
uDiff=tol+eps;

% first run
n = 1;
crit = [];
crits = tol+1;
% run until convergence or max number of iterations
while( uDiff > tol  && n < N )
    % work on other modes
    for k=1:K
        
        % recompute Hilbert mask
        HilbertMask = (sign(freqs_1*omega(n,1,k) + freqs_2*omega(n,2,k))+1);
        
        % update accumulator
        sum_uk =  sum(u_hat(:,:,:),3) - u_hat(:,:,k);
        
        % update signal spectrum
        if k==1
        u_hat(:,:,k) = ((2*rho*beta*(f_hat - sum_uk)+ u_hat(:,:,k)).*HilbertMask)./(1+2*rho*beta+ 4*rho*(Alpha(k,1)*(freqs_1 - omega(n,1,k)).^2+Alpha(k,2)*(freqs_2 - omega(n,2,k)).^2));
        else
        A = (1+2*rho*beta+ 4*rho*(Alpha(k,1)*(freqs_1 - omega(n,1,k)).^2+Alpha(k,2)*(freqs_2 - omega(n,2,k)).^2));
        B = ((2*rho*beta*(f_hat - sum_uk)+ u_hat(:,:,k)).*HilbertMask);
        param_beta = 2*(max(A(:)));
        param_fb = 1.9/(param_beta);
        un = u_hat(:,:,k);
        for l = 1:100
            %figure(3);imagesc(abs(un));axis off;colormap(gray);pause;

            un  = proj(un -  param_fb*(A.*un - B),HilbertMask);
            %un = un -  param_fb*(A.*un - B);

            un_full = fftshift(fft2(real(ifft2(ifftshift(squeeze(un))))));
            crit_fb(l) = rho*beta*norm(un_full + sum_uk - f_hat,'fro')^2 + rho*Alpha(k,2)*sum(sum((freqs_2-omega(n,2,k)).^2.*(abs(un_full).^2).*HilbertMask.^2)) + ...
            rho*Alpha(k,1)*sum(sum((freqs_1-omega(n,1,k)).^2.*(abs(un_full).^2).*HilbertMask.^2)) + 1/2*sum(sum(abs((un_full - u_hat(:,:,k)).^2)));
               
        end
  %      figure(1);
 %       semilogy(crit_fb);
%        pause;
        

        u_hat(:,:,k) = un;
        %figure(2);imagesc(abs(u_hat(:,:,k)));axis off;colormap(gray);
        end
        
        % update signal frequencies
        omega(n+1,1,k) = (omega(n,1,k) + 8*eta*Alpha(k,1)*sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2))))/(1+8*eta*Alpha(k,1)*sum(sum(abs(u_hat(:,:,k)).^2)));
        omega(n+1,2,k) = (omega(n,2,k) + 8*eta*Alpha(k,2)*sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2))))/(1+8*eta*Alpha(k,2)*sum(sum(abs(u_hat(:,:,k)).^2)));
           
        % keep omegas on same halfplane
        if omega(n+1,2,k) < 0
            omega(n+1,:,k) = -omega(n+1,:,k);
        end
        
        % recover full spectrum from analytic signal
        u_hat(:,:,k) = fftshift(fft2(real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))))));
        

    
    end
      
    % Criterion
    uDiff = eps;
    omegaDiff = eps;
    crit_tmp = beta*norm(sum(u_hat,3) - f_hat,'fro')^2;
    for k=1:K
        un_full = fftshift(fft2(real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))))));
        crit_tmp = crit_tmp  + Alpha(k,2)*sum(sum((freqs_2-omega(n+1,2,k)).^2.*(abs(un_full).^2).*HilbertMask.^2)) + ...
                    Alpha(k,1)*sum(sum((freqs_1-omega(n+1,1,k)).^2.*(abs(un_full).^2).*HilbertMask.^2))  ;
        omegaDiff = omegaDiff + sum(sum(abs(omega(n+1,:,:) - omega(n,:,:)).^2));
        uDiff = uDiff + sum(sum(1/(Hx*Hy)*(u_hat(:,:,k)-u_hat_old(:,:,k)).*conj((u_hat(:,:,k)-u_hat_old(:,:,k)))));              
                
    end
    crit(1,n) = crit_tmp;
    crit(2,n) = sum(sum(abs(sum(u_hat(:,:,:),3) - f_hat).^2));
    crits = crit(2,n);
    uDiff = abs(uDiff);
    fprintf('%d:crit =%3.2f \t crit =%3.2f\n',n,   crit(1,n), crit(2,n)); 
    n=n+1;
    u_hat_old = u_hat;
    
end


%% Signal Reconstruction

% Inverse Fourier Transform to compute (spatial) modes
u = zeros(Hy,Hx,K);
for k=1:K
    u(:,:,k) = real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))));
end;

% Should the omega-history be returned, or just the final results?
omega = omega(1:n,:,:);

end

%%
function u_hat_p = proj(u_hat,HilbertMask)
        u(:,:) = real(ifft2(ifftshift(squeeze(u_hat(:,:)))));
        %iend=floor(size(u,1)/2);
%         iend = 15;
%         iend=70;
        iend=round(size(u,1)*1/3);
        u(1:iend,:)=0;
        u_hat_p(:,:)=fftshift(fft2(u(:,:))).*HilbertMask;
end