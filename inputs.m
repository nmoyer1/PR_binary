%inputs function
% gets values for 


% Temperature
% Pressure
% Component1

function [P T x ki Tc Pc om] = inputs

% prompt_P = 'What is the  pressure in bar? >: ';
%  
prompt_T = '\n What is the flash drum  temperature in Kelvin? >: ';

prompt_frac_x1 = '\n What is the mole fraction of component one in the feed? >: '; 

%prompt_ki = '\n What is the binary interaction parameter? >:'; 

% prompt_Tc1 = '\n What is the critical temperature for component one? >';
% 
% prompt_Tc2 = '\n What is the critical temperature for component two? >:';
% 
% prompt_Pc1 = '\n What is the critical pressure for component one? >:';
% 
% prompt_Pc2 = '\n What is the critical pressure for component two? >:';
% 
% prompt_om1 = '\n What is the accentric values for component one? >:';
% 
% prompt_om2 = '\n What is the accentric values for component two? >:';

%mole fractions of feed
x1 = input(prompt_frac_x1);

x2 = 1-x1;

%Temp and pressure
%P = input(prompt_P);

P = .5; % bar
T = input(prompt_T);
 
% binary interaction parameter

%ki = input(prompt_ki);
ki = .06 % pg 441 sandler H2S and n-octane

%omega
%om1 = input(prompt_om1);
om1 = .343
%om2 = input(prompt_om2);
om2 = .592 % Perry's table 2-164 1 heptanol

%Critical Values
% Tc1 = input(prompt_Tc1);
Tc1 = 647.3 % K SAndlaer p 254
% Tc2 = input(prompt_Tc2);
Tc2 = 655 %K  NIST webbook
% Pc1 = input(prompt_Pc1);
Pc1 = 220.48 %bar Sandler
% Pc2 = input(prompt_Pc2);
Pc2 = 27 % NIST Webbook 1-octanol

x = [x1 x2];
ki = [ki 0 ; 0 ki];
om = [om1 om2];
Tc = [Tc1 Tc2];
Pc = [Pc1 Pc2];




end