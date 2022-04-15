

function [kap ai alpha bi aT xi_bi aij xi_xj_aij x_aij] = ...
    interaction(om, T, Tc ,Pc ,R, n, x, ki)



for b = 1:n
    
    %kappa vales for each species
    kap(b) = .37464 + 1.54226 * om(b) - .26992 * om(b)^2;
    
    %temp independent a values for each species
    ai(b) = .45724 *( R^2 * Tc(b)^2) / Pc(b);
    
    %temp dependent alpha values for each species
    alpha(b) = (1 + kap(b) * (1-sqrt(T / Tc(b))))^2;
    
    %temp independent b values for each species
    bi(b) = (.07780 * R * Tc(b)) / Pc(b);
    
    %temp dependent a(T) (product of ai and alpha(T)
    aT(b) = ai(b) * alpha(b); 
    
    %product of bi and mole fraction of each species for use in sum for
    %b_mix
    xi_bi(b) = x(b) * bi(b);
    
    
    
end

for i = 1:n
    for j = 1:n
        
        %aij values for a_mix calculation
        aij(i,j) = ( 1 - ki(i,j) ) * ( aT(i) * aT(j) )^.5;
        
        % product of the mole fractions and aij for the a_mix summation
        xi_xj_aij(i,j) = x(i) * x(j) * aij(i,j);
        
        %x(i) * aij for the fugacity calculation
        x_aij(i,j) = x(i) * aij(i,j);
        
        
    end
end

for i = 1:n
    for j = 1:n
        
        %aij values for a_mix calculation
        aij(i,j) = ( 1 - ki(i,j) ) * ( aT(i) * aT(j) )^.5;
        
        % product of the mole fractions and aij for the a_mix summation
        xi_xj_aij(i,j) = x(i) * x(j) * aij(i,j);
        
        %x(i) * aij for the fugacity calculation
        x_aij(i,j) = x(i) * aij(i,j);
        
        
    end
end
