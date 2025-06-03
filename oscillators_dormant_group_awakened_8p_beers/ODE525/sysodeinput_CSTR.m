% Matthew James Stephenson
% McGill ID: 261289768


function f = sysodeinput(t,x0,nvec)
%
%  This models a CSTR that may be either adiabatic nor isothermal
%  or neither.  The reaction is nA + mB --> oD .
%
%  example usage
%  sysode(m,n,xo,xf,[CA,CB,CD,CS,T])
%  sysode(1,300,0,300,[0.0,0.0,0.0,1000.0/18-0.0,500])
%
%  where
%  m = solution method
%  CA = initial concentration of A in tank [moles/liter]
%  CB = initial concentration of B in tank [moles/liter]
%  CD = initial concentration of D in tank [moles/liter]
%  CS = initial concentration of solvent in tank [moles/liter]
%  T = initial Temperature [K]
%
%  all parameters for known values 
%  and all initial guesses for unknowns 
%  are entered here
%
%
%  STEP ONE:  list the unknown variables
%
CA = x0(1);
CB = x0(2);
CD = x0(3);
CS = x0(4);
T = x0(5);
%
%  STEP TWO:  list the known parameters 
%
%  Known Flowrates [liters/time]
%  
% stream 1 (A input)
F1 = 20;	
% stream 2 (B input)
F2 = 10;
% stream C (coolant)
FC= 100;
% 
%  Known Concentrations [moles/liter]
%
% total
CT = 1000.0/18.0;
% input of A in 1
CA1 = 1.3;
% input of B in 2
CB2 = 1.7;
% input of solvent in 1
CS1 = CT - CA1;
% input of solvent in 2
CS2 = CT - CB2;
% coolant in coolant stream
CXC = CT;
%
%  Known Temperatures [K]
%
%  temp of feed 1 
T1 = 298;
%  temp of feed 2
T2 = 350;
% temp of coolant in
TCin = 273;
% temp of coolant out
TCout = 300;
if (TCout > T)
   TCout = T;
end
%
%  Reactor Volume [liters]
%
V = 1000;
%
%  Molar Heat Capacities [Joules/g-mole/K]
%
CpA = 3.0;
CpB = 5.0;
CpD = 7.0;
CpS = 4.184;
CpX = 4.184;
%
%  Volumetric Enthalpies [Joules/liter]
%
HE = T*(CpA*CA + CpB*CB + CpD*CD + CpS*CS);
H1 = T1*(CpA*CA1 + CpS*CS1);
H2 = T2*(CpB*CB2 + CpS*CS2);
HCin = TCin*CpX*CXC;
HCout = TCout*CpX*CXC;
%
%  heat input [Joules/time]
%
Qdot = FC*(HCout - HCin);
%
%  stoichiometric coefficients
%
nA = 1;
nB = 1;
nD = 2;
%
%  reaction rate constants
%
ko = 1067.1;
DHr = 109000.0;
Ea = 6203.49;
R = 8.314;
rate = (CA^nA)*(CB^nB)*ko*exp(-Ea/(R*T));

%
%  STEP THREE.  isolate the purely algebraic equations
%
CS = CT - CA - CB - CD;
FE = F1 + F2 +(nD -nA -nB)*rate/CT;
%
%  STEP FOUR .  solve for the derivatives
%
dCAdt = (F1*CA1 - FE*CA - nA*rate)/V;
dCBdt = (F2*CB2 - FE*CB - nB*rate)/V;
dCDdt = (-FE*CD + nD*rate)/V;
dCSdt = (F1*CS1 + F2*CS2 - FE*CS)/V;
num1 = F1*H1 + F2*H2 - FE*HE + DHr*rate-Qdot;
num2 = -V*T*(CpA*dCAdt + CpB*dCBdt + CpD*dCDdt + CpS*dCSdt);
den = V*(CpA*CA + CpB*CB + CpD*CD + CpS*CS);
dTdt = (num1 + num2)/den;
%
f(1) = dCAdt;
f(2) = dCBdt;
f(3) = dCDdt;
f(4) = dCSdt;
f(5) = dTdt;
%
%  STEP FIVE:  Last time through print out some values
%
if (nvec(1) == nvec(2))
   if (nvec(5) == 1)
      kstop = 1;
   else
      kstop = 4;
   end
   if (nvec(3) == kstop)
      FE
   end
end
   