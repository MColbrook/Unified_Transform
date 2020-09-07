function [a] = Gauss_Herm(n)
%compute nodes for Gauss-Lag quadrature

J=zeros(n,n);
for j=1:n
    J(j+1,j)=sqrt(j);
end

J=(J(1:n,1:n)+J(1:n,1:n)');
a=eig(J);
a=a/sqrt(2);

end


