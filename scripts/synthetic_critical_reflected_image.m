function [Usyn,Usyn_inc,Usyn_refl]=synthetic_critical_reflected_image(N,lambda,gamma,gammacm,u_inc,grid_y,piv_grid,nx,nz)

x_ph = (0:nx-1)*piv_grid/gammacm;
z_ph = (0:nz-2)*piv_grid/gammacm;

% physical parameters
nu=0.01; %viscocity [cm^2*s^-1]  
kappa=0; %cm^2*s^-1

lambda=lambda*cos(gamma);
resol=1/gammacm;
grido=resol*piv_grid;
u_max1=u_inc/(cos(2*gamma));% u_MMM is obtained from the horizontal component of the velocity
K=2*pi/lambda; % [cm^-1] wave number
phi_max=u_max1/(K);% stream function maximal value

% Spatial coordinates in real space
z_c = linspace(0,grido*grid_y,nz);
% Thierrys variables: derivate non-dimensional parameters (by order as in
% the paper). For beta=gamma
Re1=N*sin(gamma)/(nu*K^2); % Reynolds number 
Pe1=N*sin(gamma)/(kappa*K^2); % Peclet number 
a=(K^2)*sin(2*2*gamma)*phi_max/(2*N*sin(gamma)); % measure of nonlinearity
epsilon=(a*tan(2*gamma))^(1/3); 
%sigma=(gamma-gamma)/epsilon^2;
nu_6=(sin(2*gamma)^2)/(epsilon^6*Re1);
kappa_6=(sin(2*gamma)^2)/(epsilon^6*Pe1);
m2=-2/3*csc(2*gamma/2)*tan(2*gamma)*(cos(2*gamma/2)+2*(2*cos(2*gamma)-1)^(.5));
%theta=-atan(2*tan(2*gamma)/m2);
mu=cot(2*gamma/2);
rho=((nu_6+kappa_6)/(2*mu))^(.5); % some type of viscosity

% Adimensionalization of the experimental variables 

z1=K*cos(2*gamma)*z_c; % adimensional vertical coordinate

% Stretching of the experimental variables (Horizantal coordinate and vertical velocity are not stretched)
xi=tan(2*gamma)*z1/epsilon^2; % stretching vertical coordinate

% Factor of change between Thierry's and physical space velocities 
YY=(u_max1*cos(2*gamma))*tan(2*gamma)/epsilon^2;% for horizontal velocity 

% Horizontal velocity: stationary solution from Thierry's model
U_stat=zeros(1,length(z_c));
clear exp

    for j=1:length(z_c)
            if xi(j)>0
            U_stat(j)=(exp(-xi(j)/(2*rho^(2/3))))/(3^(.5)*rho^(2/3))*(sin(sqrt(3)*xi(j)/(2*rho^(2/3))+pi/3)-sqrt(3)*cos(sqrt(3)*xi(j)/(2*rho^(2/3))+pi/3));%(pi/2)*
            else
            end
    end
    
% Absolute value of the horizontal velocity
u_c=2/pi*abs(YY*U_stat); % Horizontal velocity of Thierry's solution in physical space (complete solution)



DY2D = zeros(nz,nx);
lXs = cos(gamma)*lambda/sin(2*gamma);
lZs = lXs*tan(2*gamma);
phi = 1*pi/8;
for oo=1:nx
    DY2D(:,oo) = YY*U_stat*cos(x_ph(oo)*2*pi/lXs-pi/2);
end



INC = zeros(nz,nx);

for aa = 1:length(x_ph)
    for bb = 1:length(z_c)
        INC(bb,aa) = u_inc*cos(x_ph(aa)*2*pi/lXs + z_c(bb)*2*pi/lZs-pi/2);
    end
end

Usyn = flip(DY2D+INC);
Usyn_inc = flip(INC);
Usyn_refl = flip(DY2D);

end
 
    
