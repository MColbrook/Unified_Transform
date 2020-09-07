function [a] = New_Orthog_poly_zeros(n)

load('coefficients.mat')

J=zeros(n,n);
for j=1:n
    J(j,j)=a(j)/2;
    J(j+1,j)=b(j+1);
end

J=(J(1:n,1:n)+J(1:n,1:n)');
a=eig(J);

end