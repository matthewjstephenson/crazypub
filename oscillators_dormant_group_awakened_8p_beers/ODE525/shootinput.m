% Matthew James Stephenson
% McGill ID: 261289768

function dydt = shootinput(x,y,nvec);
%dydt = zeros(o,1);
%
%  enter in equations
%
id = 2;
if (id == 1)
a = 1.0;
b = -0.1;
c = 0.01;
dydt(1) = y(2);
dydt(2) = a*x+b*y(1)+c*y(2);
elseif (id == 2)
   %
   %  homework 14, problem 4
   %  sample usage:
   % shooting(2,100,0,1,[1,0],1)
   c=[1.0,-2.0,2.0,0.0];
   dydt(1) = y(2);
   dydt(2) = c(1)*y(2) + c(2)*y(1) + c(3)*sin(x) + c(4);
elseif (id == 3)
   %
   %  homework 14, problem 4, better alternate
   %
   c=[-1.0,-2.0,2.0,0.0];
   dydt(1) = y(2);
   dydt(2) = c(1)*y(2) + c(2)*y(1) + c(3)*sin(x) + c(4);
 end