% Matthew James Stephenson
% McGill ID: 261289768


function sysode(m,n,xo,xf,yo);
%  sysode(m,n,xo,xf,yo)
%  This routine solves one non-linear first-order ordinary differential
%  equation initial value problem.
%
%  m = 1 for Euler's method
%  m = 2 for Classical Runge-Kutta 4rth order method
%  n = number of steps
%  xo = starting value of x
%  xf = ending value of x
%  o = number of first order ordinary differential equations
%  yo = initial condition at xo
%
%  The differential equation must appear in the file 'sysodeinput.m'
%  This program creates an output data file 'sysode.out'
%
%  Author:  David Keffer   Date:  October 23, 1998
%

%
%  STEP ONE.  input parameters
%
if (nargin ~= 5)
   error('sysode requires 5 input arguments');
end
nsize=size(yo);
o=max(nsize(1),nsize(2));
%
%  STEP TWO DEFINE ODE and I
%
%  yo is the initial condition
%  the ode is defined in a file called odeivpn.m
%yo = zeros(1,o);
%
%  the vector nvec contains some info that can be passed 
%  to sysodeinput (in case you need it)
%  nvec(1) = i, the iteration number
%  nvec(2) = n, the total number of iterations
%  nvec(3) = kk, the intra-iteration access number
%  nvec(4) = o, the number of first order ODEs
%  nvec(5) = m, the method used to solve
%
%  STEP THREE.  SOLVE THE ODE
%
h = (xf-xo)/n;
x = xo : h : xf;
if (m == 1) 
   y = eulerevaln(x,h,n,o,yo);
else
   y = rk4evaln(x,h,n,o,yo);
end
%
%  STEP FOUR.  PLOT THE RESULT
%
figure(1)
for i = 1:o 
  if (i==1)
     plot (x,y(:,i),'k-'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==2)
     plot (x,y(:,i),'r-'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==3)
     plot (x,y(:,i),'b-'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==4)
     plot (x,y(:,i),'g-'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==5)
     plot (x,y(:,i),'m-'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==6)
     plot (x,y(:,i),'k:'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==7)
     plot (x,y(:,i),'r:'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==8)
     plot (x,y(:,i),'b:'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==9)
     plot (x,y(:,i),'g:'), xlabel( 'x' ), ylabel ( 'y' )
  elseif (i==10)
     plot (x,y(:,i),'m:'), xlabel( 'x' ), ylabel ( 'y' )
 else
     plot (x,y(:,i),'k-'), xlabel( 'x' ), ylabel ( 'y' )
  end
  hold on
end
hold off
%
%  STEP FIVE.  WRITE THE RESULT TO sysode.out
%
op1=o+1;
fid = fopen('sysode.out','w');
if (op1 == 2) 
  fprintf(fid,'x              y(1) \n');
  fprintf(fid,'%13.7e %13.7e   \n', [x;y']);
elseif (op1 == 3)
  fprintf(fid,'x              y(1)           y(2) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e  \n', [x;y']);
elseif (op1 == 4)
  fprintf(fid,'x              y(1)           y(2)           y(3) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e %13.7e \n', [x;y']);
elseif (op1 == 5)
  fprintf(fid,'x              y(1)           y(2)           y(3)           y(4) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e %13.7e %13.7e \n', [x;y']);
elseif (op1 == 6)
  fprintf(fid,'x              y(1)           y(2)           y(3)           y(4)           y(5) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e  %13.7e %13.7e %13.7e \n', [x;y']);
elseif (op1 == 7)
  fprintf(fid,'x              y(1)           y(2)           y(3)           y(4)           y(5)           y(6) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e  \n', [x;y']);
elseif (op1 == 8)
  fprintf(fid,'x              y(1)           y(2)           y(3)           y(4)           y(5)           y(6)           y(7) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e  \n', [x;y']);
elseif (op1 == 9)
  fprintf(fid,'x              y(1)           y(2)           y(3)           y(4)           y(5)           y(6)           y(7)           y(8) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e \n', [x;y']);
elseif (op1 == 10)
  fprintf(fid,'x              y(1)           y(2)           y(3)           y(4)           y(5)           y(6)           y(7)           y(8)           y(9) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e  \n', [x;y']);
elseif (op1 == 11)
  fprintf(fid,'x              y(1)           y(2)           y(3)           y(4)           y(5)           y(6)           y(7)           y(8)           y(9)           y(10) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e  %13.7e %13.7e \n', [x;y']);
else
  fprintf(fid,'x              y(1)           y(2)           y(3)           y(4)           y(5)           y(6)           y(7)           y(8)           y(9)           y(10)           y(11) \n');
  fprintf(fid,'%13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e  \n', [x;y']);
end
fclose(fid);
%
%  If you want, interpolate to find a time at a particular value of y
%
interp = 0;
if (interp == 1)
   final_height = y(n+1)
   yf=2.54;
   time_at_h_1inch = INTERP1(y,x,yf,'linear')
end

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
   nvec(1)=i;
   xt = x(i);
   ytt = y(i,:);
   yt = ytt;
   nvec(3)=1;
   k1 = sysodeinput(xt,yt,nvec);
   xt = x(i)+h/2;
   for j = 1:o
      yt(j) = ytt(j)+h/2*k1(j);
   end
   nvec(3)=2;
   k2 = sysodeinput(xt,yt,nvec);
   for j = 1:o
      yt(j) = ytt(j)+h/2*k2(j);
   end
   nvec(3)=3;
   k3 = sysodeinput(xt,yt,nvec);
   xt = x(i)+h;
   for j = 1:o
      yt(j) = ytt(j)+h*k3(j);
   end
   nvec(3)=4;
   k4 = sysodeinput(xt,yt,nvec);
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
   dydt = sysodeinput(x(i),y(i,:),nvec);
   for j = 1:o
      y(i+1,j) = y(i,j)+h*dydt(j);
   end
end
