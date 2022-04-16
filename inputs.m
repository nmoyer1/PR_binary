%inputs function
% gets input values for 
% Temperature [K]
% octanol mole frac
% molar feed rate [mol/min]

%Assigns values for
%Pressure 0.5 bar
%kij = .06
%omega octanol & water
% Critical temp octanol & water [K]
% critical pressure octanol and water
function [P T z kij Tc Pc om feed_rate] = inputs;


prompt_T = ['\n What is the flash drum  temperature in Kelvin? \n' ...
    'ca~ 373K to 468K >: '];

prompt_frac_x1 = '\n What is the mole fraction of component one in the feed? >: '; 

prompt_feed_rate = '\n What is the feed rate in moles per minute? >: ';


%mole fractions of feed
z1 = input(prompt_frac_x1);

z2 = 1-z1;

%Temp and pressure

P = .5; % bar
T = input(prompt_T);
 
% binary interaction parameter

kij = .06; % pg 441 sandler H2S and n-octane

%omega

om1 = .343;

om2 = .592;  % Perry's table 2-164 1 heptanol

%Critical Values

Tc1 = 647.3; % K SAndlaer p 254
Tc2 = 655; %K  NIST webbook

Pc1 = 220.48; %bar Sandler
Pc2 = 27; % NIST Webbook 1-octanol


z = [z1 z2];
kij = [kij 0 ; 0 kij];
om = [om1 om2];
Tc = [Tc1 Tc2];
Pc = [Pc1 Pc2];
feed_rate = input(prompt_feed_rate);



end