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

prompt_frac_x1 = '\n What is the mole fraction of 1-octanol in the feed? >: '; 

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
C1 = 144.11083 ;
C2 = -13667.15667;
C3 = -16.8261;
C4 = 9.37e-18;
C5 = 6;
Pvap_C8OH = exp(C1 + C2/458.5 + C3*log(458.5) + C4^C5)*1e-5;
om1 = -1 -log10(Pvap_C8OH/27)  % Perry's table 2-164 1 heptanol
om2 = .344; %Sandler Water


%Critical Values
%Octanol
Tc1 = 655; %K  NIST webbook
Pc1 = 27; % NIST Webbook 1-octanol

%Water
Tc2 = 647.3; % K SAndlaer p 254
Pc2 = 220.48; %bar Sandler

z = [z1 z2];
kij = [kij 0 ; 0 kij];
om = [om1 om2];
Tc = [Tc1 Tc2];
Pc = [Pc1 Pc2];
feed_rate = input(prompt_feed_rate);



end