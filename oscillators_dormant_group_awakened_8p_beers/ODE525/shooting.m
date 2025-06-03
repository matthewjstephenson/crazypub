% Matthew James Stephenson
% McGill ID: 261289768


function shooting(m,n,xo,xf,yo,yf);
%  shooting(m,n,xo,xf,yo,yf)
%  This routine solves one non-linear first-order ordinary differential
%  equation boundary value problem using the shooting method.
%
%  m = 1 for Euler's method
%  m = 2 for Classical Runge-Kutta 4rth order method
%  n = number of steps
%  xo = starting value of x
%  xf = ending value of x
%  yo(1) = initial condition at xo
%  yo(2) = estimated initial derivative at xo
%  yf = final condition at xf
%
%  The differential equation must appear in the file 'shootinput.m'
%  This program creates an output data file 'shoot.out'
%
%
%  Author:  David Keffer   Date:  October 23, 1998
%

%
%  STEP ONE.  input parameters
%
if (nargin ~= 6)
   error('shooting requires 6 input arguments');
end
nsize=size(yo);
%  o = number of first order ordinary differential equations
o=max(nsize(1),nsize(2));
%
%  STEP TWO DEFINE ODE and I
%
%  yo is the initial condition
%  the ode is defined in a file called shootinput.m
%
%   define the tolerance and initialize the error
%
tol = 1.0e-6;
err = 1.0;
count = 0;
%
%  STEP THREE.  SOLVE THE ODE
%
h = (xf-xo)/n;
x = xo : h : xf;
while err > tol
   if (m == 1) 
      y = eulerevaln(x,h,n,o,yo);
   else
      y = rk4evaln(x,h,n,o,yo);
   end
   count = count+1;
   err = abs((yf - y(n+1,1))/yf);
   fprintf(1,'iteration %4i:  dy/dt guess: %13.7e, relative error = %13.7e  \n', count,yo(2),err);
   if (count == 1)
      yo1old = yo(2);
      yo(2) = -1.0;
      yo2old = yo(2);
      yf1old = y(n+1,1);
   else
      yf2old = y(n+1,1);
      yo(2) = yo1old + (yf - yf1old)/(yf2old - yf1old)*(yo2old - yo1old);
      yo1old = yo2old;
      yo2old = yo(2);
      yf1old = yf2old;
   end
end
%
%  STEP FOUR.  PLOT THE RESULT
%
for i = 1:o 
  if (i==1)
     plot (x,y(:,i),'k-'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==2)
     plot (x,y(:,i),'r-'), xlabel( 'x' ), ylabel ( 'y' )
  end
  hold on
end
hold off
%
%  STEP FIVE.  WRITE THE RESULT TO shoot.out
%
op1=o+1;
fid = fopen('shoot.out','w');
  fprintf(fid,'x              y(1)           y(2) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e  \n', [x;y']);
fclose(fid);



%
%  multivariate classical fourth order Runge-Kutta
%
function y = rk4evaln(x,h,n,o,yo)
%
%  specify initial values of y
%
y = zeros(n+1,o);
dydt = zeros(1,o);
k1 = zeros(1,o);
k2 = zeros(1,o);
k3 = zeros(1,o);
k4 = zeros(1,o);
yt = zeros(o,1);
ytt = zeros(1,o);
%
for j = 1:o
   y(1,j) = yo(j);
end
%  solve new values of dydt and y
nvec =[0;n;1;o;2];
for i =  1:n
   nvec=i;
   xt = x(i);
   ytt = y(i,:);
   yt = ytt;
   nvec(3)=1;
   k1 = shootinput(xt,yt,nvec);
   xt = x(i)+h/2;
   for j = 1:o
      yt(j) = ytt(j)+h/2*k1(j);
   end
   nvec(3)=2;
   k2 = shootinput(xt,yt,nvec);
   for j = 1:o
      yt(j) = ytt(j)+h/2*k2(j);
   end
   nvec(3)=3;
   k3 = shootinput(xt,yt,nvec);
   xt = x(i)+h;
   for j = 1:o
      yt(j) = ytt(j)+h*k3(j);
   end
   nvec(3)=4;
   k4 = shootinput(xt,yt,nvec);
   for j = 1:o
      y(i+1,j) = y(i,j)+h/6*(k1(j)+2*k2(j)+2*k3(j)+k4(j));
   end
end

%
%  multivariate Euler method
%
function y = eulerevaln(x,h,n,o,yo)
%  specify initial values of y
y = zeros(n+1,o);
dydt = zeros(1,o);
for j = 1:o
   y(1,j) = yo(j);
end
nvec =[0;n;1;o;1];
%  solve new values of dydt and y
for i =  1:n
   nvec(1) = i;
   dydt = shootinput(x(i),y(i,:),nvec);
   for j = 1:o
      y(i+1,j) = y(i,j)+h*dydt(j);
   end
end
