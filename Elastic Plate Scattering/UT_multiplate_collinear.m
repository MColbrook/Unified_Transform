function [Q1,Q2,Q3,q,eta,sol] = UT_multiplate_collinear(K0,theta,N,pos,ty,C1,C2,C3,Z1,Z2,x,y)
%%%% Computes the scattered field for collinear multiplate problem with different types of plates
%%%% INPUT PARAMETERS %%%%
% K0 : frequency of incident field
% theta : angle of incident field
% N : vector of number of basis functions in each segment
% pos : breakpoints for change of type of plate (inclduing gaps)
% ty : vector that describes the type of plate,
%      0=gap,1=stiff with square root basis,-1 stuff with Legendre poly basis,2=CC,3=FF,4=CF,5=FC (these effect the global relations and ty is of length 1 less than pos)
% C1 : first constant defining properties of plate
% C2 : second constant defining properties of plate
% C3 : rho*omega^2
% Z : points where we want to compute the scattered field
% x : vector along each plate where we want to compute the difference to
%     estimate the error etc.
% y : vector of points where we want to compute the plate displacement

%%%% OUTPUT %%%%
% Q : scattered field at points Z
% q : jump in scattered field at points given by cell array x
% eta : plate displacement at points given by cell array y
% sol : coefficients describing the solution

%% set parameters for collocation
Beta=K0/2;
scale=K0;%max(1/5,K0);
N(1)=2*N(1);
N(end)=2*N(end); % useful to have overcomplete basis for infinite intervals

M=4*sum(N);
H=haltonset(9);

MM=(1:M)/(M+1);
MM=MM(:);
K1=[-MM.^6];

K1=[K1;-1./K1]; % these are the collocation points
K1=[K1;exp(-1i*pi*MM)];
L=length(pos);

%% set up the linear system

A=zeros(length(K1),sum(N));

% first the left hand infinite interval
for j=1:N(1)
    A(:,j)=exp(-1i*Beta*pos(1)*(K1+1./K1)).*Besselhat(j/2,real(Beta*(K1+1./K1)),scale);
end
% now the right hand infinite interval
for j=1:N(end)
    A(:,j+sum(N(1:end-1)))=exp(-1i*Beta*pos(end)*(K1+1./K1)).*Besselhat(j/2,-real(Beta*(K1+1./K1)),scale);
end

% now the finite intervals, taking care to use correct basis functions and
% type of BC/unknown

ct=N(1)+1;
for aa=1:(L-1)
    if ty(aa)==0
        for j=1:N(aa+1)
            A(:,ct)=S_hat(j-1,-Beta*real(K1+1./K1),pos(aa),pos(aa+1));
            ct=ct+1;
        end
    elseif ty(aa)==1
        for j=1:N(aa+1)
            A(:,ct)=Beta/2*(K1-1./K1).*C_hat(j,-Beta*real(K1+1./K1),pos(aa),pos(aa+1));
            ct=ct+1;
        end
    elseif ty(aa)==-1
        for j=1:N(aa+1)
            A(:,ct)=Beta/2*(K1-1./K1).*leg_hat(j-1,-Beta*real(K1+1./K1),pos(aa),pos(aa+1));
            ct=ct+1;
        end
    else
        if ty(aa)==2
            load('spec_data_CC_double.mat','d')
        elseif ty(aa)==3
            load('spec_data_FF_double.mat','d')
        elseif ty(aa)==4
            load('spec_data_CF_double.mat','d')
        else
            load('spec_data_FC_double.mat','d')
        end
        for j=1:N(aa+1)
            A(:,ct)=(1+Beta/2*(K1-1./K1)*(-C1(aa)*(d(j)*2/(pos(aa+1)-pos(aa)))^4-C2(aa))).*Mode_hat(j,-Beta*real(K1+1./K1),pos(aa),pos(aa+1),ty(aa));
            ct=ct+1;
        end
        clear d
    end
end

% set up the forcing

b=0*K1;
for j=1:(L-1)
    alpha1=-K0*cos(theta)-Beta*(K1+1./K1);
    b=b-1i*K0*sin(theta)*(exp(1i*alpha1*pos(j+1))-exp(1i*alpha1*pos(j)))./(1i*alpha1);
end

%% now solve
sc = 1./sum(abs([A,b]),2);
II=find((sc>0)&(sc<Inf));
J=1:size(A,1);
J=setdiff(J,II);
A(J,:)=[];
b(J)=[];

sc = 1./sum(abs(A),1);
A = bsxfun(@times,sc,A);

warning('off')
sol=double((pinv(A))*(b));
sol=transpose(bsxfun(@times,sc,transpose(sol)));
warning('on')

%% now compute the whole field at the points Z
Q1=0*Z1;
F = @(x1,x2,y1,y2) -K0*besselh(1,1,K0*sqrt((x1-y1).^2+(x2-y2).^2))*1i/4.*(-y2)./sqrt((x1-y1).^2+(x2-y2).^2)+...
    K0*besselh(1,1,K0*sqrt((x1-y1).^2+(x2+y2).^2))*1i/4.*(y2)./sqrt((x1-y1).^2+(x2+y2).^2);
h2=0.0001;
F2 = @(x1,x2,y1,y2) (F(x1,x2,y1+h2,y2)-F(x1,x2,y1-h2,y2))/(2*h2);
F3 = @(x1,x2,y1,y2) (F(x1,x2,y1+h2,y2+h2)-F(x1,x2,y1+h2,y2-h2)+F(x1,x2,y1-h2,y2-h2)-F(x1,x2,y1-h2,y2+h2))/(4*h2^2);

for j=1:(L-1)
    if ty(j)>1
        [x0,w]=lgwt2(100,pos(j),pos(j+1));
        [x1,~]=lgwt2(100,-1,1);
        x0=x0(:);
        q1=0*x0;
        if ty(j)==2
            load('spec_data_CC_double.mat','d')
        elseif ty(j)==3
            load('spec_data_FF_double.mat','d')
        elseif ty(j)==4
            load('spec_data_CF_double.mat','d')
        else
            load('spec_data_FC_double.mat','d')
        end
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*Mode_val(aa,x1(:),ty(j))*(-C1(j)*(d(aa)*2/(pos(j+1)-pos(j)))^4-C2(j));
        end
        clear d
        Z=zeros(length(Z1),100);
        for aa=1:100
            Z(:,aa)=F(x0(aa),0,real(Z1(:)),imag(Z1(:)));
        end
        Q1=Q1+Z*(w(:).*q1)/2;
    elseif ty(j)==1
        [x0,w]=lgwt2(100,pos(j),pos(j+1));
        [x1,~]=lgwt2(100,-1,1);
        x0=x0(:);
        q1=0*x0;
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*C_val(aa,x1(:));
        end   
        Z=zeros(length(Z1),100);
        for aa=1:100
            Z(:,aa)=F(x0(aa),0,real(Z1(:)),imag(Z1(:)));
        end
        Q1=Q1+Z*(w(:).*q1)/2;
    elseif ty(j)==-1
        [x0,w]=lgwt2(100,pos(j),pos(j+1));
        [x1,~]=lgwt2(100,-1,1);
        x0=x0(:);
        q1=0*x0;
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*legendreP(aa-1,x1(:));
        end   
        Z=zeros(length(Z1),100);
        for aa=1:100
            Z(:,aa)=F(x0(aa),0,real(Z1(:)),imag(Z1(:)));
        end
        Q1=Q1+Z*(w(:).*q1)/2;
    end        
end

Q2=0*Z2;

for j=1:(L-1)
    if ty(j)>1
        [x0,w]=lgwt2(1000,pos(j),pos(j+1));
        [x1,~]=lgwt2(1000,-1,1);
        x0=x0(:);
        q1=0*x0;
        if ty(j)==2
            load('spec_data_CC_double.mat','d')
        elseif ty(j)==3
            load('spec_data_FF_double.mat','d')
        elseif ty(j)==4
            load('spec_data_CF_double.mat','d')
        else
            load('spec_data_FC_double.mat','d')
        end
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*Mode_val(aa,x1(:),ty(j))*(-C1(j)*(d(aa)*2/(pos(j+1)-pos(j)))^4-C2(j));
        end
        clear d
        Z=zeros(length(Z2),1000);
        for aa=1:1000
            Z(:,aa)=F2(x0(aa),0,real(Z2(:)),imag(Z2(:)));
        end
        Q2=Q2+Z*(w(:).*q1)/2;
    elseif ty(j)==1
        [x0,w]=lgwt2(1000,pos(j),pos(j+1));
        [x1,~]=lgwt2(1000,-1,1);
        x0=x0(:);
        q1=0*x0;
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*C_val(aa,x1(:));
        end   
        Z=zeros(length(Z2),1000);
        for aa=1:1000
            Z(:,aa)=F2(x0(aa),0,real(Z2(:)),imag(Z2(:)));
        end
        Q2=Q2+Z*(w(:).*q1)/2;
    elseif ty(j)==-1
        [x0,w]=lgwt2(1000,pos(j),pos(j+1));
        [x1,~]=lgwt2(1000,-1,1);
        x0=x0(:);
        q1=0*x0;
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*legendreP(aa-1,x1(:));
        end   
        Z=zeros(length(Z2),1000);
        for aa=1:1000
            Z(:,aa)=F2(x0(aa),0,real(Z2(:)),imag(Z2(:)));
        end
        Q2=Q2+Z*(w(:).*q1)/2;
    end        
end

Q3=0*Z2;

for j=1:(L-1)
    if ty(j)>1
        [x0,w]=lgwt2(1000,pos(j),pos(j+1));
        [x1,~]=lgwt2(1000,-1,1);
        x0=x0(:);
        q1=0*x0;
        if ty(j)==2
            load('spec_data_CC_double.mat','d')
        elseif ty(j)==3
            load('spec_data_FF_double.mat','d')
        elseif ty(j)==4
            load('spec_data_CF_double.mat','d')
        else
            load('spec_data_FC_double.mat','d')
        end
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*Mode_val(aa,x1(:),ty(j))*(-C1(j)*(d(aa)*2/(pos(j+1)-pos(j)))^4-C2(j));
        end
        clear d
        Z=zeros(length(Z2),1000);
        for aa=1:1000
            Z(:,aa)=F3(x0(aa),0,real(Z2(:)),imag(Z2(:)));
        end
        Q3=Q3+Z*(w(:).*q1)/2;
    elseif ty(j)==1
        [x0,w]=lgwt2(1000,pos(j),pos(j+1));
        [x1,~]=lgwt2(1000,-1,1);
        x0=x0(:);
        q1=0*x0;
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*C_val(aa,x1(:));
        end   
        Z=zeros(length(Z2),1000);
        for aa=1:1000
            Z(:,aa)=F3(x0(aa),0,real(Z2(:)),imag(Z2(:)));
        end
        Q3=Q3+Z*(w(:).*q1)/2;
    elseif ty(j)==-1
        [x0,w]=lgwt2(1000,pos(j),pos(j+1));
        [x1,~]=lgwt2(1000,-1,1);
        x0=x0(:);
        q1=0*x0;
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*legendreP(aa-1,x1(:));
        end   
        Z=zeros(length(Z2),1000);
        for aa=1:1000
            Z(:,aa)=F3(x0(aa),0,real(Z2(:)),imag(Z2(:)));
        end
        Q3=Q3+Z*(w(:).*q1)/2;
    end        
end


%% compute jump in q on the plates
q=cell(L-1,1);
for j=1:L-1
    x0=x{j};
    x0=x0(:);
    q1=0*x0;
    if ty(j)>1
        if ty(j)==2
            load('spec_data_CC_double.mat','d')
        elseif ty(j)==3
            load('spec_data_FF_double.mat','d')
        elseif ty(j)==4
            load('spec_data_CF_double.mat','d')
        else
            load('spec_data_FC_double.mat','d')
        end
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*Mode_val(aa,x0,ty(j))*(-C1(j)*(d(aa)*2/(pos(j+1)-pos(j)))^4-C2(j));
        end
        clear d
        q{j}=q1;
    elseif ty(j)==1
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*C_val(aa,x0);
        end   
       q{j}=q1;
   elseif ty(j)==-1
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*legendreP(aa-1,x0);
        end   
       q{j}=q1;
    end
end

%% compute plate displacements
eta=cell(L-1,1);
for j=1:L-1
    x0=y{j};
    x0=x0(:);
    q1=0*x0;
    if ty(j)>1
        if ty(j)==2
            load('spec_data_CC_double.mat','d')
        elseif ty(j)==3
            load('spec_data_FF_double.mat','d')
        elseif ty(j)==4
            load('spec_data_CF_double.mat','d')
        else
            load('spec_data_FC_double.mat','d')
        end
        for aa=1:N(j+1)
            q1=q1+sol(sum(N(1:j))+aa)*Mode_val(aa,x0,ty(j));
        end
        clear d
        eta{j}=q1/C3(j);
    else
        eta{j}=q1;
    end
end
sol=sol(N(1)+1:sum(N(1:end-1)));

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

function [X]=C_hat(n,K,a,b)
    K2=(b-a)/2*K;
    X=exp(1i*(a+b)/2*K).*(-n*pi./K2).*besselj(n,-K2)*(b-a)/2;
end

function [X]=leg_hat(n,K,a,b)
    K2=(b-a)/2*K;
    X=exp(1i*(a+b)/2*K).*sqrt(2*pi*K2*1i)./(K2*1i).*besseli(n+1/2,1i*K2)*(b-a)/2;
end

function [X]=S_hat(n,K,a,b)
    K2=(b-a)/2*K;
    X=exp(1i*(a+b)/2*K).*(1i)^n.*besselj(n,K2)*pi*(b-a)/2;
end

function [X]=Mode_hat(n,K,a,b,ty)
if ty==2
    load('spec_data_CC_double.mat','d','L1','L2','L3','L4','R1','R2','R3','R4')
elseif ty==3
    load('spec_data_FF_double.mat','d','L1','L2','L3','L4','R1','R2','R3','R4')
elseif ty==4
    load('spec_data_CF_double.mat','d','L1','L2','L3','L4','R1','R2','R3','R4')
else
    load('spec_data_FC_double.mat','d','L1','L2','L3','L4','R1','R2','R3','R4')
end

K2=(b-a)/2*K;
X=((exp(1i*K2)*R1(n))-(exp(-1i*K2)*L1(n))).*(1i*K2).^3-((exp(1i*K2)*R2(n))-(exp(-1i*K2)*L2(n))).*(1i*K2).^2+...
    ((exp(1i*K2)*R3(n))-(exp(-1i*K2)*L3(n))).*(1i*K2)-((exp(1i*K2)*R4(n))-(exp(-1i*K2)*L4(n)));
X=X./(K2.^4-d(n)^4);
X=(b-a)/2*exp(1i*K*(a+b)/2).*X;
end

function [X]=Mode_val(n,x,ty)
if ty==2
    load('spec_data_CC_double.mat','eta')
elseif ty==3
    load('spec_data_FF_double.mat','eta')
elseif ty==4
    load('spec_data_CF_double.mat','eta')
else
    load('spec_data_FC_double.mat','eta')
end

X=Chebfun0(250,x)*eta(1:251,n);
end

function [X]=C_val(n,x)
if mod(n,2)==1
    X=cos(n*asin(x));
else
    X=1i*sin(n*asin(x));
end
end

function [X] = Chebfun0(n,z)
    X=(zeros(length(z),n+1));
    X(:,1)=X(:,1)+1;
    X(:,2)=(z(:));
    for j=3:n+1
        X(:,j)=2*z(:).*X(:,j-1)-X(:,j-2);
    end
end
































