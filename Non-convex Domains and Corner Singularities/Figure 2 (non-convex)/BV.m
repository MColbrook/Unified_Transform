function [LD,LN] = BV(zs,ze,a,b,M)

% Create Legendre coefficients for the Dirichlet and Neumann data for the 
% test function f(z) = exp(a+b*z) along the line segment from zs to ze.

% Input parameters
%   zs,ze   Start and end points of line segment (complex)
%     a,b   Parameters defining the test function f(z) = exp(a+b*z)
%       M   Number of Legendre coefficients; use degrees up through M-1
% Output parameters
%   LD,LN   Column vectors with the first M Legendre coefficients for the
%           test function's Dirichlet and Neumann data, respectively

% Exact integral of exp(alpha*t)*P_m(t), {t,-1,1}
LI  = @(m,alpha) sqrt(2*pi*alpha)./alpha.*besseli(m+0.5,alpha);

m   = (0:M-1)';    % Column vector with the Legendre degrees to be used
LD  = (m+0.5)*exp(a+b*zs+0.5*(ze-zs)*b).*LI(m,0.5*(ze-zs)*b);
LN  = -1i*(ze-zs)*LD/abs(zs-ze);  

    