% PR flash for a binary system
% takes user input of liquid feed mole fraction, Temp [K], feed rate and
% returns a flash calculation for octanol(1) and water(2)
% for a flash drum at 0.5 bar 

%PREOS General form
% P=((R*T)/(V-b))-((a*alpha(T))/(V*(V+b)+b*(V-b)))

%Component definitions
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


% Component 1 is 1-octanol
% Component 2 is water
%general constants

% gas constant
R = 8.3145e-5; % [bar * m^3 / mol * K]

%number of components
n = 2;



% Get user inputs and assign parameters for water and octanol

[P, T, z, kij, Tc, Pc, om, feed_rate] = inputs;

% Calculate PR parameters and build interaction matrices 
[kap, ai, alpha ,bi, aT, xi_bi, aij, xi_xj_aij, x_aij] = ...
    interaction(om,T, Tc, Pc, R, n, z, kij);

% Calculate compressibility, fugacity, and fug. coeff at given T & P
[Z, k, phi_liq, phi_vap,fug_liq,fug_vap] = ...
    Z_phi_k(P, T, R,kap, ai, alpha ,bi, aT, xi_bi, aij, xi_xj_aij, x_aij, n );


flash_vec = flash(z, n, k);

