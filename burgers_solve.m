% Returns the velocity field and distance for 1D Burger's Equation
function [u, x] = burgers_solve(nt, nx, tmax, xmax, v)
% Differentials
dt = tmax/(nt-1);
dx = xmax/(nx-1);

% Initialize the computation
u = zeros(nx, nt);
x = zeros(1, nx);
ipos = x;
ineg = x;

for i = 1:nx
   x(i) = i*dx;
   ipos(i) = i+1;
   ineg(i) = i-1;
end

% Represent the periodic Boundary Conditions
ipos(nx) = 1;
ineg(1) = nx - 1;

% Represent the Initial Conditions
for i = 1:nx
   phi = exp(-(x(i)^2)/(4*v)) + exp(-(x(i)-2*pi)^2 / (4*v));
   dphi = -(0.5*x(i)/v)*exp(-(x(i)^2) / (4*v)) - (0.5*(x(i)-2*pi) / v)*exp(-(x(i)-2*pi)^2/(4*v));
   u(i,1) = -2*v*(dphi/phi) + 4;
end

% Compute the numerical solution using the discretized Burger's equation
% that was derived in my paper.
   for n = 1:nt-1
       for i = 1:nx
           u(i,n+1) = (u(i,n)-u(i,n)*(dt/dx)*(u(i,n)-u(ineg(i),n))+ ...
                      v*(dt/dx^2)*(u(ipos(i),n)-2*u(i,n)+u(ineg(i),n)));
       end
   end
end