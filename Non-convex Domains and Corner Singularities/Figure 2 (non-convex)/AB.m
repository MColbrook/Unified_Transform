function [RD,RN,SD,SN] = AB(zs,ze,k,M)

% Calculate the matrix blocks that correspond to a line segment that
% goes from point zs to point ze

%  Input parameters
% zs,ze   Start and end points of line segment (complex)
%     k  Column vector (complex), all K different k-values
%     M  Number of Legendre coefficients, i.e. degrees 0, 1, ... , M-1
%  Output parameters
%    RD  Array (K,M); part 'Regular  Dirichlet' of system matrix
%    RN  Array (K,M); part 'Regular  Neumann'   of system matrix
%    SD  Array (K,M); part 'Schwartz Dirichlet' of system matrix
%    SN  Array (K,M); part 'Schwartz Neumann'   of system matrix

% Exact integral of exp(alpha*t)*P_m(t), {t,-1,1}
LI  = @(m,alpha) sqrt(2*pi*alpha)./alpha.*besseli(m+0.5,alpha);

K = length(k); 
RD = zeros(K,M); RN = zeros(K,M);SD = zeros(K,M);SN = zeros(K,M);

for m = 0:M-1           % Loop over the degrees of Legendre polynomials
    RI = 0.5*exp(-0.5i*k*     (zs+ze))  .* LI(m,-0.5i*k*    (ze-zs));  
    RS = 0.5*exp( 0.5i*k*(conj(zs+ze))) .* LI(m, 0.5i*k*conj(ze-zs));
    RD(:,m+1) = k *    (ze-zs) .*RI;
    RN(:,m+1) =     abs(ze-zs) .*RI;    
    SD(:,m+1) = k *conj(ze-zs) .*RS;    
    SN(:,m+1) =     abs(ze-zs) .*RS;    
end  