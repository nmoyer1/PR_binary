%PR binary system

%PREOS General form
% P=((R*T)/(V-b))-((a*alpha(T))/(V*(V+b)+b*(V-b)))

%Component definitions

% Necessary inputs
% 
% Tc, T, kappa, Pc, P, omega(accentric factor) 
%  a = (.45724*R^2*Tc^2)/Pc (units of v^2/mol^2)
%  b = (.07780*R*Tc)/Pc  (units of vol/mol)
%  kappa = .37464 + 1.54226*omega - .26992*omega^2
%  alpha(T) = (1-kappa*(1-sqrt(T/Tc)))^2
%  Dimensionless form w/r/t Z (compressibility factor)
%  A = (alpha(T)*a*P)/(R^2*T^2)
%  B = (b*P)/(R*T)
%  0 = Z^3 - (1-B)*Z^2 + (A-2*B-3*B^2)*Z - (A*B - B^2 -B^3)

clc, clear all

%general comstants

% gas constant
R = 8.3145e-5; % [bar * m^3 / mol * K]
%number of components
n = 2;

[P T x ki Tc Pc om] = inputs;



[kap ai alpha bi aT xi_bi aij xi_xj_aij x_aij] = interaction(om,T, Tc, Pc, R, n, x, ki);


%PREOS combining rules 
a_mix = sum(sum(xi_xj_aij));

b_mix = sum(xi_bi);


%Calculating the compressibility factor from the a_mix and b_mix

%Coefficients
A = (a_mix * P)/(R^2 * T^2);

B = (b_mix * P)/(R * T);

%Coefficients for each degree of Z in the cubic eqn
Z3 = 1;
Z2 = -(1 - B);
Z1 = (A - 2 * B - 3 * B^2);
Z0 = -(A * B - B^2 - B^3);

% find roots of cubic polynomial 
z = [Z3 Z2 Z1 Z0];
Z = roots(z);

%delete middle value with no physical meaning
% change order to Z[liquid vapor]
Z = Z([3 1])

%calculate molar volume of liquid and vapor
Vm = (Z * R * T) / P
