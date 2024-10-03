function [M2, P2] = MP(M1, P1, theta, gamma)
% Mach number (M2) & Pressure (P2) after oblique shock
beta = obliquerelations('mach', M1, 'theta', theta, gamma); % * (180/pi) for rad to degree if it's needed
Mn1 = M1*sin(beta);
Mn2 = sqrt( (1 + ((gamma-1)/2)*Mn1^2) / (gamma*Mn1^2 - (gamma-1)/2) );
M2 = Mn2 / sin(beta-theta);
P2_P1 = 1 + ((2*gamma)/(gamma+1)) * (Mn1^2-1);
P2 = P2_P1 * P1;
