function [E1,E2] = Spec_method(N,K0)

%%%%%% INPUT PARAMETERS %%%%%%
% N - number of basis functions in each co-ordinate direction
% K0 - wavenumber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=N;
Beta=K0/2;

scale1=sqrt(2*K0);
scale2=sqrt(2*K0);

x1 = Gauss_Herm(2*N+1);
x1=x1(:)/scale1;
x1=x1(x1>0);
x1=sort(x1);

x2 = New_Orthog_poly_zeros(M);
x2=x2(:)/scale2;
x2=[0;x2];
x2=sort(x2);

fun = @(n,x) HermEval(n,x).*exp(-x.^2/2);
fun2 = @(n,x) New_Orthog_poly(n,x).*exp(-x.^2/2);

fun_d1 = @(n,x) New_Orthog_poly_d1(n,x).*exp(-x.^2/2);
fun_d2 = @(n,x) New_Orthog_poly_d2(n,x).*exp(-x.^2/2);

%% Compute differential matrix

A1=zeros(length(x1),N);
A2=A1;

for j=1:N
    A1(:,j)=fun(2*j-1,x1*scale1).*(-2*(2*j-1+1/2)*scale1^2);
    A2(:,j)=fun(2*j-1,x1*scale1);
end

B1=zeros(length(x2),M);
B2=B1;

for j=1:M
    B1(:,j)=scale2^2*fun_d2(j-1,x2*scale2)-2*x2.*scale2^3.*fun_d1(j-1,x2*scale2)-scale2^2.*fun2(j-1,x2*scale2);%-4*K0^2*x2.^2.*fun_d2(j-1,x2*scale2);
    B2(:,j)=fun2(j-1,x2*scale2);
end

A=kron(A1,B2)+kron(A2,B1);


%% Compute boundary matrix

bd=zeros(1,M);
for j=1:M
    bd(j)=scale2*fun_d1(j-1,0);
end

BD=kron(A2,bd);

%% Set up RHS and the system of equations

A=[A;BD];
b=zeros(length(x1)*length(x2)+length(x1),1);
b((length(x1)*length(x2)+1):end)=2*x1.*exp(-x1.^2/2);

sc = 1./sum(abs(A),2);
A  = bsxfun(@times,sc,A);
b  = bsxfun(@times,sc,b);
sol1 = real(A\b);

%% compute error
xd=0.01:0.01:10;
xd=xd(:);

if abs(Beta-1/4)>10^(-12)
    qd=-exp(-xd/2)/sqrt(1/4-4*Beta^2).*erfi(sqrt(1/2-2*Beta)*sqrt(xd));
else
    qd=-2*exp(-xd/2)/sqrt(2*Beta+1/2).*sqrt(xd/pi);
end

bd2=zeros(1,M);
for j=1:M
    bd2(j)=fun2(j-1,0);
end
bd3=zeros(length(xd),N);
for j=1:N
    bd3(:,j)=fun(2*j-1,sqrt(xd(:))*scale1);
end
bd4=kron(bd3,bd2);
qd_sol=bd4*sol1;

xn=0.01:0.01:1;
xn=xn(:);

qn=-1./(sqrt(pi*xn*(2*Beta+1/2))).*exp(-2*Beta*xn)+exp(xn/2).*erfc(sqrt(xn*(2*Beta+1/2)));

bd2=zeros(length(xn),M);
for j=1:M
    bd2(:,j)=fun2(j-1,sqrt(xn(:))*scale1)*sqrt(K0/2)./sqrt(xn);
end
bd3=zeros(1,N);
for j=1:N
    bd3(j)=fun(2*j-2,0)*sqrt(j-1/2)*2;
end
bd4=kron(bd3,bd2);
qn_sol=bd4*sol1;

E1=max(abs(qd-qd_sol)./max(abs(qd)));
E2=max(abs(qn-qn_sol)./max(abs(qn)));


end

