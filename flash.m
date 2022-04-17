
function flash_vec = flash(z, n ,k)
    
syms x1 x2 L

eqn1 = x1 + x2 == 1;
eqn2 = x1*(L*(1-k(1))+k(1)) == z(1);
eqn3 = x2*(L*(1-k(2))+k(2)) == z(2);

soln = solve([eqn1, eqn2, eqn3], [x1, x2, L])
double(soln.x1)
double(soln.x2)
double(soln.L)
liq_split = double(soln.L(1))
vap_split = 1-liq_split
split_sum = liq_split + vap_split;

C8OH_L = double(soln.x1(1))
H2O_L = double(soln.x2(1))
%liq_sum = (acetone_L+acetonitrile_L);

C8OH_V = C8OH_L * k(1)
H2O_V = H2O_L * k(2)
%vap_sum = (C8OH_V+H2O_V);

flash_vec = [C8OH_L H2O_L C8OH_V H2O_V liq_split vap_split ]

end

