% Matthew James Stephenson
% McGill ID: 261289768

function dydt = sysodeinput(x,y,nvec);
%
%  This is the CO2 problem
%
%  sample inputs 
%
%   sysode(2,1000,0,10,[1.0;  0.0; 0.0; 0.0; 0.0; 0.0])
%   sysode(2,1000,0,10,[1.0;  1.0; 1.0; 0.0; 0.0; 0.0])
%   sysode(2,1000,0,10,[0.0;  0.0; 0.0; 1.0; 0.0; 0.0])
%   sysode(2,1000,0,10,[0.0;  0.0; 0.0; 1.0; 1.0; 1.0])
%   sysode(2,1000,0,10,[0.0;  0.0; 0.0; 1.0; 0.0; -1.0])
%
   m1 = 16;
   m2 = 12;
   m3 = 16;
   k12 = 48;
   k13 = 48;
   A = [ 	0				0					0					1					0					0
      		0				0					0					0					1					0
				0				0					0					0					0					1
            -k12/m1		k12/m1			0					0					0					0
            k12/m2		(-k12-k13)/m2	k13/m2			0					0					0
            0				k13/m3			-k13/m3			0					0					0];
         
        
   dydt = A*y';