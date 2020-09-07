function [I] = sing_int(rho,alpha,m)
    if m==0
        I =1./(-rho).^((alpha+1)).*(gamma(alpha+1)-igamma(alpha+1,-2*rho));
    else
        V=0;
        for j=0:m
            V=V+nchoosek(m,j)*(log(2))^(m-j)*factorial(j)*(-1)^(j)/((alpha+1)^(j+1))*hypergeom(repmat(alpha+1,1,j+1), repmat(alpha+2,1,j+1), 2*rho);
        end
        I=2^(alpha+1)*V;
    end
end