close all
clear
clc
%% Solving Q5
% Initial Conditions
gamma = 1.4;
M(1) = 3;
theta(2) = 9 * (pi/180);
theta(1) = 5 * (pi/180);
theta(3) = theta(1);
P(1) = 50; % (kPa)

% 2 & 3 Area
[M(2), P(2)] = MP(M(1), P(1), theta(2), gamma);
[M(3), P(3)] = MP(M(1), P(1), theta(3), gamma);

% Solving Shock System
phi_g = 2 * (pi/180); % Guessed Phi
d = 10 * (pi/180); % Calculation domain (-/+)10 degree
phi = phi_g;
[M(4), P(4)] = MP(M(2), P(2), theta(2) + phi, gamma);
[M(5), P(5)] = MP(M(3), P(3), theta(1) - phi, gamma);
% Iterating
tol = 1e-6;
while abs(P(4)-P(5)) > tol
    if P(4) == P(5)
        break
    elseif P(4) < P(5)
            phi = phi + d;
            [M(4), P(4)] = MP(M(2), P(2), theta(2) + phi, gamma);
            [M(5), P(5)] = MP(M(3), P(3), theta(1) - phi, gamma);
    else
            phi = phi - d;
            [M(4), P(4)] = MP(M(2), P(2), theta(2) + phi, gamma);
            [M(5), P(5)] = MP(M(3), P(3), theta(1) - phi, gamma);
    end
    d = d/2;
end
P4 = P(4)
P5 = P(5)
Phi = phi * (180/pi)