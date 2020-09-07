function [X] = laguerreL2(Arg1,Arg2,Arg3)

if nargin <3
    x=Arg2;
    n=Arg1;
    a=1;
    b=1-x;
    if n<0
        X=0*x;
    elseif n==0
        X=a;
    elseif n==1
        X=b;
    else
        for j=2:n
            c=(2*j-1-x).*b-(j-1)*a;
            c=c/(j);
            a=b;
            b=c;
        end
        X=b;
    end
else
    x=Arg3;
    alpha=Arg2;
    n=Arg1;
    a=1;
    b=1-x+alpha;
    if n<0
        X=0*x;
    elseif n==0
        X=a;
    elseif n==1
        X=b;
    else
        for j=2:n
            c=(2*j-1-x+alpha).*b-(j-1+alpha)*a;
            c=c/(j);
            a=b;
            b=c;
        end
        X=b;
    end
end

end

