function [u_analytical, x] = analytical_solution(nt, nx, tmax, xmax, nu)

% Discretization of increments
dt = tmax/(nt-1);
dx = xmax/(nx-1);

% Initialise data structures
u_analytical = zeros(nx, nt);
x = zeros(1, nx);
t = zeros(1, nt);

% Computes the displacement
for i = 1:nx
   x(i) = i*dx;
end

% Computes the Analytical Solution to test accuracy
for n = 1:nt
   t = n*dt;

   for i = 1:nx
       phi = exp(-(x(i)-4*t)^2/(4*nu*(t+1))) + exp(-(x(i)-4*t-2*pi)^2/(4*nu*(t+1)));

       dphi = (-0.5*(x(i)-4*t)/(nu*(t+1))*exp(-(x(i)-4*t)^2/(4*nu*(t+1)))...
           -0.5*(x(i)-4*t-2*pi)/(nu*(t+1))*exp(-(x(i)-4*t-2*pi)^2/(4*nu*(t+1))));

       u_analytical(i,n) = -2*nu*(dphi/phi) + 4;
   end

end
end