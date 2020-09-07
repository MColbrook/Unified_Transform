function [X] = M_hankel(parity,q,N,m,xi)
%%% computes M Hankel function

[~,B] = Eig_problem(parity,q,N);
B=B(:,m);

v1=sqrt(q)*exp(-xi);
v2=sqrt(q)*exp(xi);

X=0*xi;
if parity==0
    for j=1:N
        X=X+(-1)^j*B(j)*(besselj(j-1,v1).*besselh(j+1,1,v2)- besselj(j+1,v1).*besselh(j-1,1,v2)   );
    end
else
    for j=1:N
        X=X-(-1)^j*B(j)*(besselj(j-1,v1).*besselh(j,1,v2)- besselj(j,v1).*besselh(j-1,1,v2)   );
    end
end

X=X/M_hankel_d(parity,q,N,m,0);  % normalisation
end

