% Matthew James Stephenson
% McGill ID: 261289768

%% prepare a random disordered metric and dilaton factor 
N=201;
n = [9 9]; %%% number of points in each element
N0 = [40,40]; %%% number of elements in each direction

NGlobal = [N,N];
%[dx,dy,d2x,d2y,xx,yy] = mkSpectralOperators({N,'fourier'},{N,'fourier'});
xx = linspace(-pi, pi, N+1);
xx = xx(1:N);
yy = xx;
%del = {dx,dy};

[x,y] = ndgrid(xx,yy); 

nrand = floor(N/20);
epsrand = 0.05;   %disorder strength in metric
zrand = 0.07;     %disorder strength in dilaton factor

gxx = zeros(N);
gxy0 = zeros(N);
gyy = zeros(N);
Z = zeros(N);

for ii=1:nrand
    for jj=1:nrand
        gxx = gxx + (2*rand()-1)*epsrand*cos(2*pi*rand()+ii*x).*cos(2*pi*rand()+jj*y);
        gxy0 = gxy0 + (2*rand()-1)*epsrand*cos(2*pi*rand()+ii*x).*cos(2*pi*rand()+jj*y);
        gyy = gyy + (2*rand()-1)*epsrand*cos(2*pi*rand()+ii*x).*cos(2*pi*rand()+jj*y);
        Z = Z + (2*rand()-1)*zrand*cos(2*pi*rand()+ii*x).*cos(2*pi*rand()+jj*y);
    end
end

gxx = exp(gxx);
gyy = exp(gyy);
gxy = 1-exp(-gxy0.^2);
gxy = sign(gxy0).*(gxx.*gyy.*gxy).^(1/2);
Z = exp(Z);

detg = sqrt(gxx.*gyy-gxy.^2);

ggxx = Z.*gyy./detg;
ggyy = Z.*gxx./detg;
ggxy = -1*Z.*gxy./detg;

%% prepare to use different code...

dims = 2;
%parpool(4)

C2 = cell(dims,dims);
C1 = cell(dims,dims);
C2{1,1} = ggxx;
C2{2,2} = ggyy;
C2{1,2} = ggxy;
C2{2,1} = ggxy;
C1{1}   = derx(ggxx,N,N) + dery(ggxy,N,N);
C1{2}   = dery(ggyy,N,N) + derx(ggxy,N,N);
C0      = 0*x;

Sx = -C1{1};
Sy = -C1{2};


alphax = SolveLinearSystem(n,N0,NGlobal,C2,C1,C0,Sx);
alphay = SolveLinearSystem(n,N0,NGlobal,C2,C1,C0,Sy);

%surf(xx,yy,alphax),view(0,90), xlim([-pi pi]), ylim([-pi pi]), colorbar, shading interp, lighting phong, shg

sigmaxx = mean(mean(ggxx.*(1+derx(alphax,N,N)) + ggxy.*dery(alphax,N,N) ))
sigmaxy = mean(mean(ggxy.*(1+dery(alphay,N,N)) + ggxx.*derx(alphay,N,N) ))
sigmayx = mean(mean(ggxy.*(1+derx(alphax,N,N)) + ggyy.*dery(alphax,N,N) ))
sigmayy = mean(mean(ggyy.*(1+dery(alphay,N,N)) + ggxy.*derx(alphay,N,N) ))
detsigma = sqrt(sigmaxx*sigmayy-sigmaxy*sigmayx)

%% set up Poisson equation on metric
%note -- slow part of algorithm is computing L2 and sources!   Can do this for sure because we are awesome
%by using fast fourier transform for derivatives acting on ggxx though...

% L2 = Dx*(ggxx*Dx+ggxy*Dy)+Dy*(ggxy*Dx+ggyy*Dy);
% 
% 
% L2 = L2(2:N^2,2:N^2);
% 
% u = ones(N*N,1);
% 
% sourcex = -Dx*ggxx*u-Dy*ggxy*u;
% sourcey = -Dx*ggxy*u-Dy*ggyy*u;
% 
% sourcex = sourcex(2:N^2);
% sourcey = sourcey(2:N^2);
% 
% alphax = L2\sourcex;
% alphay = L2\sourcey;
% 
% alphax = [0; alphax];
% alphay = [0; alphay];
% 
% 
% sigmaxx = mean((ggxx*Dx+ggxy*Dy)*alphax + ggxx*u)
% sigmaxy = mean((ggxy*Dx+ggyy*Dy)*alphax + ggxy*u)
% sigmayx = mean((ggxx*Dx+ggxy*Dy)*alphay + ggxy*u)
% sigmayy = mean((ggxy*Dx+ggyy*Dy)*alphay + ggyy*u)
% detsigma = sqrt(sigmaxx*sigmayy-sigmaxy*sigmayx)