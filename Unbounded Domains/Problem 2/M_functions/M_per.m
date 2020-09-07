function [X] = M_per(parity,q,N,m,eta)
%%% computes M periodic functions (odd)


[~,B] = Eig_problem(parity,q,N);
B=B(:,m);

X=0*eta;
if parity==0
    for j=1:N
        X=X+B(j)*sin(2*j*eta);
    end
else
    for j=1:N
        X=X+B(j)*sin((2*j-1)*eta);
    end
end

end
