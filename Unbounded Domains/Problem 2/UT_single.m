function [q1] = UT_single(K0,theta,N1)
%%% compute approximation via UT method

%% set parameters
a=1/2;
Beta=K0/2;
scale=K0;
N3=N1;
N2=2*N1;
N1=2*N1;

M1=3*N1;
M2=M1;
M3=M1;

H=haltonset(9);
K1=[-H(2:M1+1,1);
    -H(2:M2+1,2).^2*0.1];
K1=[K1;exp(-1i*pi*H(2:M3+1,5))];
K1=[K1;-1./K1];

%%

A=zeros(length(K1),N1+N2+N3);
A2=zeros(length(K1),N1);

for j=1:N1
    A2(:,j)    = Besselhat(j/2,real( Beta*(K1+1./K1)),scale);
end
A(:,1:N1)    = A2(:,1:N1).*repmat(exp(a*1i*Beta*(K1+1./K1)),1,N1);
A(:,(N1+1):(N1+N2))=conj(A2(:,1:N1)).*repmat(exp((a-1)*1i*Beta*(K1+1./K1)),1,N2);

for j=1:N3
    A(:,j+N1+N2)    = Beta/4*(K1-1./K1).*Chebhat(j,real(-Beta/2*(K1+1./K1))).*exp(1i*(a-1/2)*Beta*(K1+1./K1));
end


theta=-theta;
alpha1=K0*cos(theta)-Beta*(K1+1./K1);
b=1i*K0*sin(theta)*(exp(1i*(1-a)*alpha1)-exp(-1i*a*alpha1))./(1i*alpha1);

%% now implicitly include zero end-point BCs

A(length(K1)+(1:10),(N1+N2+2):2:(N1+N2+N3))=repmat(sin(pi/2*(2:2:(N3))),10,1);
b(length(K1)+(1:10))=zeros(10,1);

sc = 1./sum(abs([A,b]),2);
I=find((sc>0)&(sc<Inf));
J=1:size(A,1);
J=setdiff(J,I);
A(J,:)=[];
b(J)=[];

sc = 1./sum(abs(A),1);
A = bsxfun(@times,sc,A);

sol=double((pinv(A))*(b));
sol=transpose(bsxfun(@times,sc,transpose(sol)));

%% reconstruct the solution

sol1=sol((N1+N2+1):end);

t=-1:0.01:1;
q1=0*t;

for j=1:N3
    q1=q1+sol1(j)*Chebval(j,t);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X] = Besselhat(mu,z,scale)
X=(z>0).*( (z>scale).*scale^mu.*exp(mu/2*pi*1i)./(mu*(z+sqrt(z.^2-scale^2)).^mu)...
    +(z<scale).*exp(1i*mu*asin(z/scale))./mu)...
    +(z<0).*conj(( ((-z)>scale).*scale^mu.*exp(mu/2*pi*1i)./(mu*((-z)+sqrt((-z).^2-scale^2)).^mu)...
    +((-z)<scale).*exp(1i*mu*asin((-z)/scale))./mu));
X=X/scale;
end

function [X] = Chebhat(n,z)
X=-n*pi*besselj(n,-z)./z;
end

function [X] = Chebval(n,z)
    if mod(n,2)==0
        X=1i*sin(n*asin(z));
    else
        X=cos(n*asin(z));
    end
end