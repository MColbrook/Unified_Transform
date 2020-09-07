function [E1,E2] = UT_method(Nd1,Nd2,Nn1,Nn2,K0)

%%%%%% INPUT PARAMETERS %%%%%%
% Nd1 - number of smooth basis functions for Dirichlet unknown
% Nd2 - number of non-smooth basis functions for Dirichlet unknown
% Nn1 - number of smooth basis functions for Neumann unknown
% Nn2 - number of non-smooth basis functions for Dirichlet unknown
% K0 - wavenumber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=Nd1+Nd2+Nn1+Nn2;
Beta=K0/2;
r2=1.015;
r1=round(K0*N*3);

% parameters for generalised Laguerre polynomials
alpha1=1/2;
alpha2=-1/2;

% There are eight rescale parameters. For each choice of basis function
% there is one for the exponential and one for the algebraic part.

scale1=1; % smooth d exp
scale2=1; % smooth d alg
scale3=1; % nonsmooth d exp
scale4=1; % nonsmooth d alg

scale5=2*K0; % smooth n exp
scale6=2*K0; % smooth n alg
scale7=2*K0; % nonsmooth n exp
scale8=2*K0; % nonsmooth n alg

H=haltonset(1);  
K=-(H(1:r1)+10*eps);
K=K(:)*r2; % these are the collocation points

%% construct and solve the linear system

L=length(K);
K2=Beta*(K-1./K);

A=zeros(L,Nn1);
B=zeros(L,Nd1);
A2=zeros(L,Nn2);
B2=zeros(L,Nd2);

for j=1:Nd1
    n=j-1;
    B(:,j)=(1/scale2)*0.5*Beta*(K+1./K).*((scale1/scale2+2i*K2/scale2-2)./(scale1/scale2+2i*K2/scale2)).^n*(2)./(scale1/scale2+2i*K2/scale2);
end
for j=1:Nd2
    n=j-1;
    B2(:,j)=(1/scale4)*0.5*Beta*(K+1./K).*(scale3/(scale4*2)+K2*1i/scale4).^(-3/2)/beta(n+1,1/2)*sqrt(pi).*(1-1./(scale3/(scale4*2)+K2*1i/scale4)).^n;
end

for j=1:Nn1
    n=j-1;
    A(:,j)=(1/scale6)*((2i*K2/scale6+2-scale5/scale6)./(2i*K2/scale6-scale5/scale6)).^n*(2)./(scale5/scale6-2i*K2/scale6);
end
for j=1:Nn2
    n=j-1;
    A2(:,j)=(1/scale8)*(scale7/(scale8*2)-K2*1i/scale8).^(-1/2)*beta(n+1/2,1/2)/sqrt(pi).*(1-1./(scale7/(scale8*2)-K2*1i/scale8)).^n;
end

b=-2./(1+2i*K2);

A=[([A, B, A2, B2]);conj([A, B, A2, B2])];
b=[(b);conj(b)];
sc = 1./sum(abs(A),2);

% Scale the system - make each row of the A matrix have 1-norm equal to one
A  = bsxfun(@times,sc,A);
b  = bsxfun(@times,sc,b);
sc = 1./sum(abs(A),1);
A  = bsxfun(@times,sc,A);


warning('off','all')
sol = real(A\b);
sol=transpose(bsxfun(@times,sc,transpose(sol)));
warning('on','all')

%% compute error
x=0.01:0.01:10;qd_sol=x*0;
for j=1:Nd2
    qd_sol=qd_sol+sol(Nd1+Nn1+Nn2+j)*laguerreL2(j-1,alpha1,x*scale4).*exp(-x*scale3/2).*sqrt(x*scale4);
end
for j=1:Nd1
    qd_sol=qd_sol+sol(Nn1+j)*laguerreL2(j-1,x*scale2).*exp(-x*scale1/2);
end
if abs(Beta-1/4)>10^(-12)
    qd=2*exp(-x/2)/sqrt(1/4-4*Beta^2).*erfi(sqrt(1/2-2*Beta)*sqrt(x));
else
    qd=4*exp(-x/2)/sqrt(2*Beta+1/2).*sqrt(x/pi);
end
E1=max(abs(real(qd-qd_sol))./(max(abs(real(qd)))));

x=0.01:0.01:1;
qn_sol=x*0;
for j=1:Nn2
    qn_sol=qn_sol+sol(Nd1+Nn1+j)*laguerreL2(j-1,alpha2,x*scale8).*exp(-x*scale7/2)./sqrt(x*scale8);
end
for j=1:Nn1
    qn_sol=qn_sol+sol(j)*laguerreL2(j-1,x*scale6).*exp(-x*scale5/2);
end
qn=-1./(sqrt(pi*x*(2*Beta+1/2))).*exp(-2*Beta*x)+exp(x/2).*erfc(sqrt(x*(2*Beta+1/2)));
E2=max(abs(real(qn-qn_sol))./(max(abs(real(qn)))));
end

