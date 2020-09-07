function [q] = M_exp(K0,theta,t,N1,N2,N3)
%%% computes expansion in M functions

%%% N1 is number of M functions
%%% N2 is truncation for e-value problem
%%% N3 is truncation for M_per functions

Q=K0^2/16;

%%% compute the coefficients for expansion
[~,B0] = Eig_problem(0,Q,N2);
[~,B1] = Eig_problem(1,Q,N2);
S=zeros(1,N1);
for j=1:N1
    if mod(j,2)==0
        S(j)=(2*M_hankel(0,Q,N1,round(j/2),0))*2i*tan(theta)*sum((1i).^(2*(1:N2)-1).*(2:2:2*N2).*transpose(B0(:,round(j/2))).*besselj(2:2:2*N2,cos(theta)*K0/2));
    else
        S(j)=(2*M_hankel(1,Q,N1,round((j+1)/2),0))*2i*tan(theta)*sum((1i).^(2*(1:N2)-2).*(1:2:2*N2).*transpose(B1(:,round((j+1)/2))).*besselj(1:2:2*N2,cos(theta)*K0/2));
    end
end

%%% compute expansion
q=0*t;
for j=1:N1
    q=q+S(j)*M_per(mod(j,2),Q,N3,floor((j+1.0001)/2),acos(t));
end

end

