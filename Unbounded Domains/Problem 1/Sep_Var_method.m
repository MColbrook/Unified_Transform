function [E1,E2] = Sep_Var_method(N,K0)

%%%%%% INPUT PARAMETERS %%%%%%
% N - number of basis functions in each co-ordinate direction
% K0 - wavenumber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Beta=K0/2;
[x1] = Gauss_Herm(2*N+1);
x1=x1(:);
x1=x1(x1>0);
x1=sort(x1);

L=length(x1);
A=zeros(L,N);

fun = @(n,x) HermEval(n,x).*exp(-x.^2/2);

%%

for j=1:N
    A(:,j)=fun(2*j-1,x1);
end

b=-1/(K0*sqrt(pi))*x1.*exp(-x1.^2/(4*K0));

sc = 1./sum(abs(A),2);
A  = bsxfun(@times,sc,A);
b  = bsxfun(@times,sc,b);

sol1 = real(A\b);
%% compute error
x=0.01:0.01:10;
qd_sol=x*0;
for j=1:N
    qd_sol=qd_sol+beta(j,1/2)*fun(2*j-1,sqrt(x*2*K0))*sol1(j)/2;
end
if abs(Beta-1/4)>10^(-12)
    qd=-exp(-x/2)/sqrt(1/4-4*Beta^2).*erfi(sqrt(1/2-2*Beta)*sqrt(x));
else
    qd=-2*exp(-x/2)/sqrt(2*Beta+1/2).*sqrt(x/pi);
end
E1=max(abs(real(qd-qd_sol))./(max(abs(real(qd)))));

x=0.01:0.01:1;
qn_sol=x*0;
for j=1:N
    qn_sol=qn_sol+sol1(j)*fun(2*j-2,0)*sqrt(j-1/2)*sqrt(2*K0)*DCyliner(j*2,2*sqrt(K0)*sqrt(x))./(2*sqrt(x));
end
qn=-1./(sqrt(pi*x*(2*Beta+1/2))).*exp(-2*Beta*x)+exp(x/2).*erfc(sqrt(x*(2*Beta+1/2)));
E2=max(abs(real(qn-qn_sol))./(max(abs(real(qn)))));

end

