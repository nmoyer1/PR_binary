function [Z, k, phi_liq, phi_vap,fug_liq,fug_vap] = ...
    Z_phi_k(P, T, R,kap, ai, alpha ,bi, aT, xi_bi, aij, xi_xj_aij, x_aij, n );

%Peng-Robinson Z and phi calculations
%takes P, R, xi_xj_aij, xi_bi
%returns Vm (molar volume), Z ( compress. factor) phi_liq, phi_vap, k1, k2 

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
Vm = (Z * R * T) / P;

liq_Vm = Vm(1);
vap_Vm = Vm(2);


%Sum of mole fractions and aij for the phi calculation
x_aij_fug = sum(x_aij);

%Fugacity coefficient calculation

for l = 1:n
    phi_liq(l) = exp(((bi(l) / b_mix) * (Z(1) - 1)) - (log(Z(1) - B)) - ...
        ((A /(2 * sqrt(2) * B)) * ((2 * x_aij_fug(l) / a_mix) - ...
        (bi(l) / b_mix)) * log((Z(1) + 2.414 * B) / (Z(1) - 0.414 * B))));

     
     phi_vap(l) = exp((bi(l) / b_mix) * (Z(2) - 1) - log(Z(2) - B) - (a_mix /...
         (2 * sqrt(2) * b_mix * R * T)) * (((2*x_aij_fug(l)) / a_mix) - (bi(l) / b_mix))...
         * log((Z(2) + (1 + sqrt(2) * B)) / (Z(2) + (1 - sqrt(2) * B))));
end

%display phi and fugacity for debugging
phi_liq;

phi_vap ;

fug_liq = phi_liq * P;

fug_vap = phi_vap * P;


%  k values using phi/phi
k1 = (phi_liq(1) / phi_vap(1));

k2 = (phi_liq(2) / phi_vap(2)) ;
k = [k1 k2];
end