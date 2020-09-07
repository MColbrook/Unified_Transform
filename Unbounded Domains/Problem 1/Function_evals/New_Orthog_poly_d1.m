function [X] = New_Orthog_poly_d1(n,x)

load('coefficients.mat')

ac=a;
bc=b;

a=0*x;
b=1+0*x;

a=a/sqrt((sqrt(pi)/2));
b=b/sqrt(0.161018670952500863350502);

if n<0
    X=0*x;
elseif n==0
    X=a;
elseif n==1
    X=b;
else
    for j=2:n
        c=((x-ac(j)).*b-bc(j)*a)/bc(j+1)+New_Orthog_poly(j-1,x)/bc(j+1);
        a=b;
        b=c;
    end
    X=b;
end

end