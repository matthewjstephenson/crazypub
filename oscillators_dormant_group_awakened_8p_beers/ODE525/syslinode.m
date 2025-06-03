% Matthew James Stephenson
% McGill ID: 261289768


function syslinode(A,yo,to,tf)
%
%  This Matlab routine will solve a system of LINEAR first-order ordinary differential equations.
%
%  syslinode(A,yo,to,tf)
%  A is the matrix of constant coefficents of the equation
%  dy/dt = A*t
%  yo is the vector of initial conditions
%  to is the value of the independent variable, t, at the initial conditions
%  tf is teh final value of interest of the independent variable, t
%
%  Author:  David Keffer, University of Tennessee, Knoxville
%  Date:  May 18, 1998
%
%  syslinode runs a sample 3x3 problem
%  syslinode('demo1') runs a first order rxn eqlbm problem
%  syslinode('demo2') runs the CO2 vibrational mode analysis problem
%  syslinode('demo3') runs a sample 2x2 problem
%

%
%  STEP -1.  CHECk AND MAKE SURE THE INPUT ARGUMENTS ARE OK
%
scrid = 1;
tol = 1.0e-6;
tol2 = 1.0e-15;
tol3 = 1.0e-6;
%  set iprint to 1 if you want to print intermediate values.
%  any other value will not print.
iprint = 0;
%  set normalize to 1 if you want the eigenvectors normalized so that their components sum to unity.
%  set normalize to 1 if you want the eigenvectors and eigenrows normalized so that they form orthonormal basis sets.
%  any other value will normalize the sum of the square of the eigenvectors to unity.
normalize = 0;
%
if (nargin ==1 & strcmp(A,'demo1') == 1) 
   [A,yo,to,tf,normalize] = get_demo_1(scrid);
elseif (nargin ==1 & strcmp(A,'demo2') == 1)
   global mass
   [A,yo,to,tf,normalize] = get_demo_2(scrid);
elseif (nargin ==1 & strcmp(A,'demo3') == 1) 
   [A,yo,to,tf,normalize] = get_demo_3(scrid);
elseif ( (nargin > 0 & nargin < 4) | (nargin > 4) )
   error('syslinode requires 4 input arguments');
elseif (nargin == 0)
   [A,yo,to,tf,normalize] = get_demo_0(scrid);
end
%
nvec = size(A);
n = nvec(1);
if (nvec(1)  ~= nvec(2))
   error('syslinode requires a square matrix to be input.');
end
nvecy = size(yo);
ny = max(nvec);
if (ny  ~= nvec(1))
   error('syslinode requires the initial condition vector to be the same length as the input matrix.');
end
detA = det(A);
rankA = rank(A);
if (iprint == 1)
   fprintf (scrid,'The starting time is %9.2e and the ending time is %9.2e \n', to, tf);
   A
   yo
	fprintf (scrid,'The determinant of A is %9.2e and the rank of A is %2i \n', detA, rankA);
end   
%
%  STEP ONE.  DETERMINE THE EIGENVALUES AND EIGENVECTORS OF THE MATRIX A
%
[wcol,lambdac] = eig(A);
for j = 1:1:n
   lam(j) = lambdac(j,j);
end
%
%  normalize the vectors if they are not normalized to sum to unity
%  (the computer normalizes them so that the sum of their squares is unity)
%
if (normalize == 1)
      fprintf(scrid,'normalizing eigenvectors\n');
   for j = 1:1:n
      wsum(j) = sum(wcol(:,j));
   end
   for j = 1:1:n
      if (abs(wsum(j)) > 1.0e-6)
      	for i = 1:1:n
         	wcol(i,j) = wcol(i,j)/wsum(j);
         end
      end
   end
end
%
%  STEP TWO.  DETERMINE THE ADJOINT OF THE MATRIX A
%
Astar = A';
%
%  STEP THREE.  DETERMINE THE EIGENVALUES AND EIGENROWS OF THE MATRIX A
%  STEP THREE.  OR DETERMINE THE EIGENVECTORS OF THE ADJOINT MATRIX A
%  (WHERE THE EIGENVALUES ARE THE SAME).
%
[wrow,lambdar]= eig(Astar);
for j = 1:1:n
   lamr(j) = lambdar(j,j);
end
%
%  normalize the vectors if they are not normalized to sum to unity
%  (the computer normalizes them so that the sum of their squares is unity)
%
if (normalize == 1)
   fprintf(scrid,'normalizing eigenrows\n');
   for j = 1:1:n
      wsum(j) = sum(wrow(:,j));
   end
   for j = 1:1:n
      if (abs(wsum(j)) > 1.0e-6)
  	    	for i = 1:1:n
         	wrow(i,j) = wrow(i,j)/wsum(j);
         end
      end
   end
end
%
%  rearrange the eigenrows if the computer did not return them in the same order as the eigenvectors
%
lamrtemp = lamr;
wrowtemp = wrow;
if (sum(abs(lamrtemp - lam)) > tol) 
   for i = 1:1:n
      for j = 1:1:n
			if ( abs(lam(i) - lamrtemp(j)) < tol )
            wrow(:,i) = wrowtemp(:,j);
            lamr(i) = lamrtemp(j);
         end
      end
   end
end
%
% STEP SIX.  DETERMINE IF YOUR EIGENVALUES ARE PURELY COMPLEX CONJUGATE PAIRS.
%
imaginary = 0;
for j = 1:1:n
   if (abs(imag(lam(j))) > tol)  
      imaginary = 1;
   end
end
%
%  STEP SEVEN.  AVOID COMPUTER ROUND-OFF AND TRUNCATION ERROR WITH COMPLEX EIGENVECTORS, VALUES, and ROWS 
%
if (imaginary == 1)
   [wrow,wcol,lam,lamr] = cleanup_imaginary(n, wrow, wcol, lam, lamr, tol3);
end   
%
%  do some testing to see if the code is working
%
for i = 1:1:n
   for j = 1:1:n
      prod(i,j) = dotpc(wrow(:,i),wcol(:,j),n);
   end
end   
if (normalize == 2)
   fprintf(scrid,'normalizing for an orthonormal basis set\n ');
   for j = 1:1:n
		denominator = dotpc( wrow(:,j), wcol(:,j),n);
   	wrow(:,j) = wrow(:,j)/denominator;
   end
end
for i = 1:1:n
   for j = 1:1:n
      prod(i,j) = dotpc(wrow(:,i),wcol(:,j),n);
   end
end   
%prod
%wrow'*wcol
%wrow'*A*wcol
%
%  STEP FOUR.  CALCULATE THE CONSTANTS WHICH SATISFY THE INITIAL CONDITIONS
%
%k = zeros(n,1);
for j = 1:1:n
   denominator = dotpc( wrow(:,j), wcol(:,j),n);
   numerator = dotpc( wrow(:,j), yo,n);
   if (imaginary == 2)
     [denominator, numerator] = cleanup_imaginary_k(n, denominator, numerator, tol3)
  	end
   if (abs(denominator) < 1.0e-12)
      k(j) = 0.0;
   else
      k(j) = numerator/denominator;
   end
end

if (imaginary == 1)
   %
   %  one of the k's doesn't exist.
   %  in that case, choose k's so that IC's are met.
   %
   yp0 = zeros(n,1);
   j = 1;
   for i = 1:2:n-3
      preexpfacr = real(k(i) * wcol(j,i));
      yp0(j) = yp0(j) +  2*preexpfacr;
   end
   k(5) = ( yo(j) - yp0(j) )/(2*wcol(j,5));
%   k(5) = k(5) + k(5)*sqrt(-1);
   k(6) = k(5);
end
%
%
%
if (imaginary == 1)
   %
   %  for second order problems, you can insert a linear term
   %
   klin = zeros(n,1);
   kcon = zeros(n,1);
   numerator = 0;
   denominator = 0;
   for j = n/2+1:1:n
      jj = j - n/2;
      numerator = numerator + mass(jj)*yo(j);
      denominator = denominator + mass(jj);
   end
   klinp = numerator/denominator;
	for j = 1:1:n/2
      klin(j) = klinp;
      kcon(j+n/2)=klinp;
   end
end
%
%  STEP FIVE.  WRITE THE SOLUTION
%
if (iprint == 1)
   k
   wcol
   lam
   wrow
   lamr
end
for j = 1:1:n
   if (abs(lam(j)) < 1.0e-15 )
      lam(j) = 0.0;
   end
end
%printout(scrid,n,k,wcol,lam,imaginary)

%
%  STEP SIX.  PLOT THE SOLUTION
%
npoints = 1000;
dt = (tf-to)/npoints;
for i = 1:1:npoints+1
   tp(i) = (i-1)*dt +to;
end
yp = zeros(n,npoints+1);
if (imaginary == 0)
   for j = 1:1:n
      for i = 1:1:npoints+1
         for ii = 1:1:n
            yp(j,i) = yp(j,i)+ k(ii) * wcol(j,ii) * exp( lam(ii) *( tp(i) - to));
         end
   	end
   end
   np = n;
else
   %
   %  create a matrix of prefactors and exponential args
   %
   nte1 = (1-to)*(npoints/tf) + 1 ;
   nte0 = (0-to)*(npoints/tf) + 1 ;
   for j = 1:1:n
      for i = 1:1:n
         pefr(j,i) = 2*real(k(i) * wcol(j,i));
         pefi(j,i) = 2*(imag(k(i) * wcol(j,i)));
         ef(j,i) = (imag(lam(i)));
      end
   end
   for j = 1:1:n
      for m = 1:1:npoints+1
        		% fprintf (scrid, 'j = %2i, m = %2i, ii = %2i, preexpfac = %9.2e, expfac = %9.2e \n ',j,m,ii,preexpfac, expfac); 
%        		yp(j,m) = yp(j,m) +  preexpfacr  * cos( expfac * ( tp(m) - to ) ) - preexpfaci  * sin( expfac * ( tp(m) - to ) ) + klin(j)*( tp(m)-to) ;
         for i = 1:2:n-1
				yp(j,m) = yp(j,m) +  pefr(j,i) * cos( ef(j,i) * ( tp(m) - to ) ) - pefi(j,i)  * sin( ef(j,i) * ( tp(m) - to ) )  ;
         end
%         yp01(j) = - yp(j,nte1) - kcon(j);
         yp(j,m) = yp(j,m) + klin(j)*( tp(m)-to)  + kcon(j);
      end
      yp00(j) =  yo(j) - (yp(j,nte0) - klin(j)*( tp(nte0)-to)  - kcon(j))  ;        
   end
   np = n;
end
%yp01
%klin
%yp00
%kcon
yp0 = yp(:,1);
iplotextra = 0;
if (iplotextra == 1)
	xbond = 0.0;
	yp(1,:) = yp(1,:) - xbond;
	yp(3,:) = yp(3,:) + xbond;
   np = np+1;
   for m = 1:1:npoints+1
      yp(np,m) = (16*yp(1,m) + 12*yp(2,m) + 16*yp(3,m))/60;   
   end
end
figure(2)
plotout(np,tp,yp)





















%
%  This old function will plot out multiple data sets in different colors
%       1 = black   solid
%       2 = red     solid
%       3 = blue    solid
%       4 = green   solid
%       5 = magenta solid
%       6 = black   dotted
%       7 = red     dotted
%       8 = blue    dotted
%       9 = green   dotted
%      10 = magenta dotted
%     >10 = black   solid
%
function plotout(n,tp,yp)
for j = 1:1:n 
  if (j==1)
     plot (tp,yp(j,:),'k-'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==2)
     plot (tp,yp(j,:),'r-'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==3)
     plot (tp,yp(j,:),'b-'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==4)
     plot (tp,yp(j,:),'g-'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==5)
     plot (tp,yp(j,:),'m-'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==6)
     plot (tp,yp(j,:),'k:'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==7)
     plot (tp,yp(j,:),'r:'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==8)
     plot (tp,yp(j,:),'b:'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==9)
     plot (tp,yp(j,:),'g:'), xlabel( 't' ), ylabel ( 'y' )
  elseif (j==10)
     plot (tp,yp(j,:),'m:'), xlabel( 't' ), ylabel ( 'y' )
 else
     plot (tp,yp(j,:),'k-'), xlabel( 't' ), ylabel ( 'y' )
  end
  hold on
end
hold off



function printout(scrid,n,k,wcol,lam, imaginary)
if (imaginary ==0)
	if (n == 1) 
   	fprintf(scrid,'x(j,t) = k(1) * wcol(j,1) * exp( lam(1) *(t-to)) \n ');
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * %9.2e * exp( %9.2e *(t-to))  \n ', j, k(1),wcol(j,1),lam(1));
   	end   
   	fprintf (scrid, ' OR \n ');
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * exp( %9.2e *(t-to)) \n ', j, k(1)*wcol(j,1),lam(1));
   	end   
	elseif (n == 2)
   	fprintf(scrid,'x(j,t) = k(1) * wcol(j,1) * exp( lam(1) *(t-to)) + k(2) * wcol(j,2) * exp( lam(2) *(t-to))  \n ');
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * %9.2e * exp( %9.2e *(t-to)) + %9.2e * %9.2e * exp( %9.2e *(t-to)) \n ', j, k(1),wcol(j,1),lam(1), k(2),wcol(j,2),lam(2));
   	end   
   	fprintf (scrid, ' OR \n ');
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * exp( %9.2e *(t-to)) + %9.2e * exp( %9.2e *(t-to))  \n ', j, k(1)*wcol(j,1),lam(1), k(2)*wcol(j,2),lam(2));
   	end   
	elseif (n == 3)
   	fprintf(scrid,'x(j,t) = k(1) * wcol(j,1) * exp( lam(1) *(t-to)) + k(2) * wcol(j,2) * exp( lam(2) *(t-to)) + k(3) * wcol(j,3) * exp( lam(3) *(t-to)) \n ');
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * %9.2e * exp( %9.2e *(t-to)) + %9.2e * %9.2e * exp( %9.2e *(t-to)) + %9.2e * %9.2e * exp( %9.2e *(t-to))  \n ', j, k(1),wcol(j,1),lam(1), k(2),wcol(j,2),lam(2), k(3),wcol(j,3),lam(3));
   	end   
   	fprintf (scrid, ' OR \n ');
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * exp( %9.2e *(t-to)) + %9.2e * exp( %9.2e *(t-to)) +  %9.2e * exp( %9.2e *(t-to))  \n ', j, k(1)*wcol(j,1),lam(1), k(2)*wcol(j,2),lam(2), k(3)*wcol(j,3),lam(3));
   	end  
	elseif (n == 4)
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * exp( %9.2e *(t-to)) + %9.2e * exp( %9.2e *(t-to)) +  %9.2e * exp( %9.2e *(t-to)) +  %9.2e * exp( %9.2e *(t-to)) \n ', j, k(1)*wcol(j,1),lam(1), k(2)*wcol(j,2),lam(2), k(3)*wcol(j,3),lam(3), k(4)*wcol(j,4),lam(4));
   	end  
	elseif (n == 5)
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * exp( %9.2e *(t-to)) + %9.2e * exp( %9.2e *(t-to)) +  %9.2e * exp( %9.2e *(t-to)) +  %9.2e * exp( %9.2e *(t-to))+  %9.2e * exp( %9.2e *(t-to)) \n ', j, k(1)*wcol(j,1),lam(1), k(2)*wcol(j,2),lam(2), k(3)*wcol(j,3),lam(3), k(4)*wcol(j,4),lam(4), k(5)*wcol(j,5),lam(5));
   	end  
	elseif (n == 6)
   	for j = 1:1:n
      	fprintf(scrid,'x(%2i,t) = %9.2e * exp( %9.2e *(t-to)) + %9.2e * exp( %9.2e *(t-to)) +  %9.2e * exp( %9.2e *(t-to)) +  %9.2e * exp( %9.2e *(t-to))+  %9.2e * exp( %9.2e *(t-to))+  %9.2e * exp( %9.2e *(t-to)) \n ', j, k(1)*wcol(j,1),lam(1), k(2)*wcol(j,2),lam(2), k(3)*wcol(j,3),lam(3), k(4)*wcol(j,4),lam(4), k(5)*wcol(j,5),lam(5), k(6)*wcol(j,6),lam(6));
   	end  
	end
else
   if (n == 6)
      fprintf(scrid,'x(j,t) = k(1) * wcol(j,1) * 2cos( lam(1) *(t-to)) + k(2) * wcol(j,2) * 2cos( lam(2) *(t-to)) + k(3) * wcol(j,3) * 2cos( lam(3) *(t-to)) \n ');   
      for j = 1:1:n/2
         preexpfac1 = real(k(1) * wcol(j,1));
         expfac1 = abs(imag(lam(1))) ;
         preexpfac3 = real(k(3) * wcol(j,3));
         expfac3 = abs(imag(lam(3))) ;
         preexpfac5 = real(k(5) * wcol(j,5));
         expfac5 = abs(imag(lam(5))) ;
      	fprintf(scrid,'x(%2i,t) = %9.2e *cos( %9.2e *(t-to)) + %9.2e *cos( %9.2e *(t-to)) +  %9.2e *cos( %9.2e *(t-to)) \n ', j, preexpfac1,expfac1,  preexpfac3,expfac3,  preexpfac5,expfac5  );
   	end
   end
end






function    [wrow,wcol,lam,lamr] = cleanup_imaginary(n, wrow, wcol, lam, lamr, tol3)
tol4 = 1.0e-7;
%
%  clean up eigenrows
%
  for m = 1:1:n
    	for j = 1:1:n
   		if (abs(real(wrow(m,j))) < tol3*abs(imag(wrow(m,j))) )
      		realpart = 0.0;
   		else
         	realpart = real(wrow(m,j));
      	end
   		if (abs(imag(wrow(m,j))) < tol3*abs(real(wrow(m,j))) )
      		imagpart = 0.0;
   		else
         	imagpart = imag(wrow(m,j));
         end
         wrow(m,j) = realpart + imagpart*sqrt(-1);
         if ( abs(real(wrow(m,j))) < tol4 & abs(imag(wrow(m,j))) < tol4 )
            wrow(m,j) = 0.0;
         end
      end  
   end
%
%  clean up eigenvectors
%
   for m = 1:1:n
    	for j = 1:1:n
   		if (abs(real(wcol(m,j))) < tol3*abs(imag(wcol(m,j))) )
      		realpart = 0.0;
   		else
         	realpart = real(wcol(m,j));
      	end
   		if (abs(imag(wcol(m,j))) < tol3*abs(real(wcol(m,j))) )
      		imagpart = 0.0;
   		else
         	imagpart = imag(wcol(m,j));
      	end
         wcol(m,j) = realpart + imagpart*sqrt(-1);
         if ( abs(real(wcol(m,j))) < tol4 & abs(imag(wcol(m,j))) < tol4 )
            wcol(m,j) = 0.0;
         end
      end  
   end
%
%  clean up eigenvaluess
%
   for m = 1:1:n
   		if (abs(real(lam(m))) < tol3*abs(imag(lam(m))) )
      		realpart = 0.0;
   		else
         	realpart = real(lam(m));
      	end
   		if (abs(imag(lam(m))) < tol3*abs(real(lam(m))) )
      		imagpart = 0.0;
   		else
         	imagpart = imag(lam(m));
      	end
         lam(m) = realpart + imagpart*sqrt(-1);
         if ( abs(real(lam(m))) < tol4 & abs(imag(lam(m))) < tol4 )
            lam(m) = 0.0;
         end
         lamr(m) = lam(m);  
   end

   
   
   function   [A,yo,to,tf,normalize] = get_demo_1(scrid)
   fprintf(scrid,'demo problem #1 is the first-order rxn eqlbm problem \n');
   k12 = 0.50;
   k21 = 0.25;
   k13 = 0.20;
   k31 = 0.05;
   k23 = 0.30;
   k32 = 0.15;
   A = [ 	(-k13-k12)	k21				k31
      k12					(-k21-k23)		k32
      k13					k23				(-k31-k32)];
   yo = [1.0/3.0;  1.0/3.0; 1.0/3.0];
   to = 0.0;
   tf = 10.0;
   normalize = 1;
   
   function   [A,yo,to,tf,normalize] = get_demo_2(scrid)
   fprintf(scrid,'demo problem #2 is the CO2 equations of motion problem \n');
   global mass
   mass(1) = 16;
   mass(2) = 12;
   mass(3) = 16;
   k12 = 48;
   k13 = 48;
   A = [ 	0					0							0					1					0					0
      		0					0							0					0					1					0
				0					0							0					0					0					1
            -k12/mass(1)	k12/mass(1)				0					0					0					0
            k12/mass(2)		(-k12-k13)/mass(2)   k13/mass(2)		0					0					0
            0					k13/mass(3)				-k13/mass(3)	0					0					0];
   yo = [1.0;  0.0; 0.0; 1.0; 0.0; 0.0];
   to = 0.0;
   tf = 10;
   normalize = 0;
   
   function   [A,yo,to,tf,normalize] = get_demo_3(scrid)
   fprintf(scrid,'demo problem #3 is just a test problem\n');
   k12 = 0.50;
   k21 = 0.25;
   A = [ 	-2*k12	k21
      		k12 	-k21];
   yo = [1.0/2.0;  1.0/2.0];
   to = 5.0;
   tf = 10.0;
   normalize = 0;
   
   function   [A,yo,to,tf,normalize] = get_demo_0(scrid)
   fprintf(scrid,'Okay, we will evaluate a sample problem \n');
   A = [ 2 -2 3; 1 1 1; 1 3 -1 ];
   yo = [1.0;  0.0; 0.5];
   to = 0.0;
   tf = 1.0;
   normalize = 0;
   
   
   
   
         function   [denominator, numerator] = cleanup_imaginary_k(n, denominator, numerator, tol3)
         real(denominator)
         imag(denominator)
         real(numerator)
         imag(numerator)
         if ( abs(real(denominator)) < tol3)
            realden = 0.0;
         else
            realden = real(denominator);
         end
         if ( abs(imag(denominator)) < tol3)
            imagden = 0.0;
         else
            imagden = imag(denominator);
         end
         denominator = realden + imagden*sqrt(-1);
         if ( abs(real(numerator)) < tol3)
            realnum = 0.0;
         else
            realnum = real(numerator);
         end
         if ( abs(imag(numerator)) < tol3)
            imagnum = 0.0;
         else
            imagnum = imag(numerator);
         end
         numerator = realnum + imagnum*sqrt(-1);
         
         
      function [answer] = dotpc(a, b, n)
      answer = 0.0; 
      for i = 1:1:n
   		answer = answer + a(i)*b(i);
		end   
