function [X] = M_hankel_d(parity,q,N,m,xi)
%%% compute derivative of M Hankel function pre-normalisation

[~,B] = Eig_problem(parity,q,N);
B=B(:,m);

v1=sqrt(q)*exp(-xi);
v2=sqrt(q)*exp(xi);

X=0*xi;
if parity==0
    for j=1:N
        X=X+(-1)^j*B(j)*(  -v1.*(besselj(j-2,v1)-besselj(j,v1))/2.*besselh(j+1,1,v2)+v2.*besselj(j-1,v1).*(besselh(j,1,v2)-besselh(j+2,1,v2))/2 ... 
                           +v1.*(besselj(j,v1)-besselj(j+2,v1))/2.*besselh(j-1,1,v2)-v2.*besselj(j+1,v1).*(besselh(j-2,1,v2)-besselh(j,1,v2))/2);
    end
else
    for j=1:N
        X=X-(-1)^j*B(j)*(  -v1.*(besselj(j-2,v1)-besselj(j,v1))/2.*besselh(j,1,v2)+v2.*besselj(j-1,v1).*(besselh(j-1,1,v2)-besselh(j+1,1,v2))/2 ...
                          +v1.*(besselj(j-1,v1)-besselj(j+1,v1))/2.*besselh(j-1,1,v2)-v2.*besselj(j,v1).*(besselh(j-2,1,v2)-besselh(j,1,v2))/2   );
    end
end
end