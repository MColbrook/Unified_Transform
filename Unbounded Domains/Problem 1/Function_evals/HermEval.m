function [X] = HermEval(n,x)
%Evaluate Hermite functions recursively (without exponential weight)

a=1/(pi^0.25);
b=sqrt(2)*x.*a;

if n<0
    X=0*x;
elseif n==0
    X=a;
elseif n==1
    X=b;
else
    for j=2:n
        c=-sqrt((j-1)/j)*a+sqrt(2/j)*x.*b;
        a=b;
        b=c;
    end
    X=b;
end

end

