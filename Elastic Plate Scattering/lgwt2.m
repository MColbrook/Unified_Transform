function  [x,w]=lgwt2(N,a,b)
[x,w]=legpts(N); % Chebfun routine, computes Legendre points and Gauss-Legendre quadrature weights.
% Linear map from[-1,1] to [a,b]
x=a+(b-a)*(x+1)/2;    
% Compute the weights
w=(b-a)/2*w;