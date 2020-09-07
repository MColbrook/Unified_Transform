function [b,B] = Eig_problem(parity,q,N)
%%% solves the eigenvalue problem needed to build M functions

A=zeros(N);
if parity==0
    for j=1:N-1
        A(j,j+1)=q;
        A(j+1,j)=q;
    end
    for j=1:N
        A(j,j)=(2*j)^2;
    end
else
    for j=1:N-1
        A(j,j+1)=q;
        A(j+1,j)=q;
    end
    for j=1:N
        A(j,j)=(2*j-1)^2;
    end
    A(1,1)=A(1,1)-q;
end
    
[V,D] = eig(A);
b=diag(D);
B=V;
end

