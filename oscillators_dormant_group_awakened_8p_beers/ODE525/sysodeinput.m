% Matthew James Stephenson
% McGill ID: 261289768

function dydt = sysodeinput(x,y,nvec);
%
%  sample input file for sysode
%
%  example usage:
%  sysode(2,100,0,1,[-1,1])
%
%  enter in equations
%
a = 1.0;
b = -0.1;
c = 0.011;
dydt(1) = y(2);
dydt(2) = a*x+b*y(1)/(c*y(2))-y(2);