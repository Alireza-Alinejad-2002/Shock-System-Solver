% Function to solve the mach, beta, theta oblique shock relations for
% supersonic flow incident on a wedge.
%
% Currently works for solving weak shock angles
%
% Mach - the mach number of the flow (dimensionless)
% Beta - the shock angle (radians)
% Theta - the wedge half-angle (in radians).  In some cases called the turning angle of
%         the flow, but this code is not intended to solve Prandtle-Meyer expansion
%         problems.
%
% See, e.g., http://en.wikipedia.org/wiki/Oblique_shock for a description
% of the geometric relationship between beta and theta.
%
% Given a mach number and beta, theta is solved using 
% 
%  tan(theta) = 2cot(beta)*(M^2sin^2(beta) - 1)/(M^2(gamma + cos(2beta) + 1)
%
% Given a theta and a mach number, the above equation is solved numerically
% using Newton's method. 
%
% Given a theta and a beta, the above equation is solved for M using
% Newton's method.
%
% The two input vectors must be of equal length.  If they are not, the
% shorter input vector will be lengthened and the last value of the shorter
% vector will be repeated to expand the vector.
function retval = obliquerelations(label1, in1, label2, in2, gamma)
retval = -1;
% Expand the shortest vector to be the same length as the longest vector.
% Hold the last value of the shorter vector.
if length(in1) < length(in2)
    newin1 = zeros(size(in2));
    newin1(1:length(in1)) = in1;
    newin1(length(in1):length(newin1)) = ones(size(in2))*in1(length(in1));
    in1 = newin1;
elseif length(in1) > length(in2)
    newin2 = zeros(size(in1));
    newin2(1:length(in2)) = in2;
    newin2(length(in2):length(newin2)) = ones(size(in1))*in2(length(in2));
    in2 = newin2;
end
if length(in1) ~= length(in2)
    error('Inputs must be of the same length');
else
    % Given a mach and a beta, solve for theta directly
    if strcmpi(label1, 'mach') && strcmpi(label2, 'beta')
        retval = atan2(2*cot(in2).*(in1.^2.*(sin(in2)).^2.-1), in1.^2.*(gamma + cos(2*in2))+2);
    elseif strcmpi(label2, 'mach') && strcmpi(label1, 'beta')
        retval = atan2(2*cot(in1).*(in2.^2.*(sin(in1)).^2.-1), in1.^2.*(gamma + cos(2*in1))+2);
        
    % Given a theta and a beta, solve for M using Newton's method
    elseif strcmpi(label1, 'beta') && strcmpi(label2, 'theta')
        b = in1;
        theta = in2;
        f = @(M) (2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2)-tan(theta);
        fp = @(M) -1*(4*M.*cot(b).*(M.^2.*sin(b).^2-1).*(cos(2*b)+gamma))./((M.^2.*(cos(2*b)+gamma)+2).^2)...                
                + (4*M.*cos(b).*sin(b))./(M.^2.*(cos(2*b)+gamma) +2);
            
        xold = ones(size(in1))*1.5;
        xnew = xold - f(xold)./fp(xold);
        while abs(xold - xnew) > 10^-8*ones(size(in1))
            
            xold = xnew;
            xnew = xold - f(xold)./fp(xold);
        end
        retval = xnew;
    elseif strcmpi(label1, 'theta') && strcmpi(label2, 'beta')
        b = in2;
        theta = in1;
        f = @(M) (2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2)-tan(theta);
        fp = @(M) -1*(4*M.*cot(b).*(M.^2.*sin(b).^2-1).*(cos(2*b)+gamma))./((M.^2.*(cos(2*b)+gamma)+2).^2)...                
                + (4*M.*cos(b).*sin(b))./(M.^2.*(cos(2*b)+gamma) +2);
             
        xold = ones(size(in1));
        xnew = xold - f(xold)./fp(xold);
        while abs(xold - xnew) > 10^-8*ones(size(in1))
            
            xold = xnew;
            xnew = xold - f(xold)./fp(xold);
        end
        retval = xnew;
        
    % Given a mach and a theta, solve for beta using Newton's method
    % (return NaN for out of bounds values)
    elseif strcmpi(label2, 'mach') && strcmpi(label1, 'theta')
        M = in2;
        theta = in1;
        f = @(b) (2*cot(b).*(M.^2.*(sin(b)).^2-1))/(M.^2.*(gamma + cos(2*b))+2)-tan(theta);
        fp = @(b) (4*M.^2.*sin(2*b).*cot(b).*(M.^2.*sin(b).^2-1))./((M.^2.*(cos(2*b)+gamma)+2).^2)...                
                + (4*M.^2.*cos(b).^2 - 2*csc(b).^2.*(M.^2.*(sin(b)).^2-1))./(M.^2.*(cos(2*b)+gamma) +2);
        overmax = ones(size(theta));
        for i = 1:length(in1)
           M = in1(i);
           f1 = @(b) -1*(2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2);
           x = fminbnd(f1,0,pi/2);
           if(theta(i) >= atan(-1*f1(x)))
               overmax(i) = NaN;
           end
        end  
        xold = ones(size(in1))*0.1;
        xnew = xold - f(xold)/fp(xold);
        while abs(xold - xnew) > 10^-8*ones(size(in1))
            
            xold = xnew;
            xnew = xold - f(xold)/fp(xold);
        end
        retval = xnew;
        retval = retval.*overmax;
    elseif strcmpi(label1, 'mach') && strcmpi(label2, 'theta')
        M = in1;
        theta = in2;
        f = @(b) (2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2)-tan(theta);
        fp = @(b) (4*M.^2.*sin(2*b).*cot(b).*(M.^2.*sin(b).^2-1))./((M.^2.*(cos(2*b)+gamma)+2).^2)...                
                + (4*M.^2.*cos(b).^2 - 2*csc(b).^2.*(M.^2.*(sin(b)).^2-1))./(M.^2.*(cos(2*b)+gamma) +2);
        overmax = ones(size(theta));
        for i = 1:length(in1)
           M = in1(i);
           f1 = @(b) -1*(2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2);
           x = fminbnd(f1,0,pi/2);
           if(theta(i) >= atan(-1*f1(x)))
               overmax(i) = NaN;
           end
        end
        xold = ones(size(in1))*0.1;
        xnew = xold - f(xold)./fp(xold);
        while abs(xold - xnew) > 10^-8*ones(size(in1))
            xold = xnew;
            xnew = xold - f(xold)./fp(xold);
        end
        retval = xnew;
        retval = retval.*overmax;
    end
end