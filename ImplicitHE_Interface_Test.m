function [tmax,xmax] = ImplicitHE_Interface_test(dx,dt)

%Written by James Rekow
%Based on the code ImplicitHE written by Dr. Vrushali Bokil

% solve the heat equation using implicit scheme (equation for v)
% solve Vt=Vxx for x in (0,1), t>0
%  V(0,t)=V(1,t)=0, V(x,0)=f(x)= 1 (when x<0.5), 0 (otherwise)

%End time
T=0.1; 

%all the inner grid points
m=round(1/dx)-1; N=round(T/dt); r=dt/(dx*dx);
x=[dx:dx:1-dx];

%b will be the initial data
b=zeros(size(x))';

%Identifies which grid point corresponds to (0.5+)
k = 1;

%initial condition:
for i=1:m,

  if x(i)<0.5,
      k = k + 1;
      if x(i)>0.2,
          if x(i)<0.4,
            b(i) = exp(-(1/(1-0.01*(x(i)-0.3)^2)));
          end
      else
          b(i) = 0;
      end
  
  else 
    b(i)=0;
  end,
end,

%This will be the x value at which u attains its maximum
xmax = x(k);

%maximum value of b(k)
yv = b(k);
yu = 10*yv;
tmax = 0;

%initial data
c = b;

%solution u to the equivalent interface problem
u = b;

for j = k:length(x),
    u(j) = 10*b(j);
end
   

% The implicit scheme gives us the system of linear equations
%             (I-kA)U^{n+1} = U^n
% since the coefficient matrix B = I-kA is tri-diagonal with constant 
% sup/sub/ diagonal elements, we can use three scalar values to 
% present A, namely, alpha, beta and gamma. 
alpha=1+2*r; beta=-r; gamma=-r;

% Build the matrix B
e = ones(m,1);
B = alpha*diag(e,0)+beta*diag(e(2:m),-1)+gamma*diag(e(1:m-1),1);

%time
t = 0;

%Vectors used to create solution surface
tvec = zeros(1,N);
uvec = zeros(N,length(x));

for nT=1:N,
    
% Solve the system by Gaussian elimination
b = B\b;

%
u = b;
for j = k:length(x),
    u(j) = 10*b(j);
end

for r=1:length(x),
    uvec(nT,r) = u(r);
end

%update max value at x = 0.5
if b(k) > yv;
    yv = b(k);
    yu = 10*yv;
    tmax = t;
end
  

figure(1)

  %plot solution to transformed problem and initial data
subplot(2,1,1)
plot(x,c,'--',x,b,'--o',x(k),yv,'*');
axis([0,1,0,5]);
h=text(0.2,0.05,sprintf('dx=%g, dt=%g, r=%g,t=%g',dx,dt,r,t)); 
set(h,'FontSize',12),

% add title, labels, and legends to the above plot
title('Solution to Transformed Problem')
xlabel('x')
ylabel('v')


  %plot solution to interface problem and initial data
subplot(2,1,2)
plot(x,c,'--',x,u,'--o',x(k),yu,'*');
axis([0,1,0,5]);
h=text(0.2,0.05,sprintf('dx=%g, dt=%g, r=%g,t=%g',dx,dt,r,t)); 
set(h,'FontSize',12),

% add title, labels, and legends to the above plot
title('Solution to Interface Problem')
xlabel('x')
ylabel('u')

t = t + dt;
tvec(nT) = t;

pause(0.0003)
  
end

[X,Y] = meshgrid(x,tvec);

figure(2)

surf(X,Y,uvec)
view(-45,45)
hold on;
plot3(xmax,tmax,yu,'sk','markerfacecolor',[1,0,0]);

% add title, labels, and legends to the above plot
title('Solution Surface u(x,t) to Interface Problem')
xlabel('x')
ylabel('t')
zlabel('u')
