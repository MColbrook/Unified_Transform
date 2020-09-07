function [ cjl,AA] = Lap_sing( N,col,S )
% N number of Leg polys
% col the collocation parameters
% S the number of singular functions
% now leave out the non singular part along relevant edges

vertices=[0,1,1+2i,1i];
n=4;
for j=1:3
    hj(j)=0.5*(vertices(j+1)-vertices(j));
    mj(j)=0.5*(vertices(j+1)+vertices(j));
end
hj(4)=0.5*(vertices(1)-vertices(4));
mj(4)=0.5*(vertices(1)+vertices(4));


if strcmp(col.type,'ray')
    R1=col.R1;
    R2=col.R2;
    M=col.M;
    lvec=linspace(R1,R2,M); 

    k=[];

    for j=3:4
        fvec=-lvec(1:end)/hj(j);
        k=[k(:);fvec(:)];
    end
    for j=1:2
        fvec=-(lvec(2:2:end))/hj(j);
        k=[k(:);fvec(:)];
    end
    
    R=col.R;
    K=col.K;
    H   = haltonset(2);  % 2-D node set
    kr  = 2*H(1:2*K,:)-1;
    kt  = R*complex(kr(:,1),kr(:,2)); % Convert to complex
    kt   = kt(find(abs(kt)<R,K));
    k=[k(:);kt(:);conj(kt(:))];
elseif strcmp(col.type,'halton')

    R=col.R;
    K=col.K;
    H   = haltonset(2);  % 2-D node set
    kr  = 2*H(1:2*K,:)-1;
    kt  = R*complex(kr(:,1),kr(:,2)); % Convert to complex
    k   = kt(find(abs(kt)<R,K));
end
        
% Create the matrices that form the system matrix
z=vertices;
[~,RN1,~,~] = AB(z(1),z(2),k,N);          % Side 1
[RD2,~,~,~] = AB(z(2),z(3),k,N);          % Side 2

% compute contribution from singular functions


Ncoef=-2/3:(-4/3):-1000;
Ncoef=Ncoef(1:S);
RSIN=cell(1,S);
RSID=RSIN;
I=1:length(k);

digitsOld=digits(70);
for j=1:S
    RSIN{j}=zeros(length(k),1);
    RSID{j}=zeros(length(k),1);  
    RSIN{j}(I)=Ncoef(j)*sing_int(vpa(-1i*hj(4)*k(I)),vpa(-Ncoef(j)-1),0).*exp(vpa(1i*hj(4)*k(I)-1i*mj(4)*k(I)))*(abs(hj(4)))^(-Ncoef(j)-1)*abs(hj(4));
    RSID{j}(I)=sing_int(vpa(1i*hj(3)*k(I)),vpa(-Ncoef(j)),0).*exp(vpa(-1i*hj(3)*k(I)-1i*mj(3)*k(I)))*(abs(hj(3)))^(-Ncoef(j)).*(k(I)*hj(3))*sin(3*pi/4*(-Ncoef(j)));
end
digits(digitsOld);

%%
BB = -[ RN1, RD2];

SING=[];
for j=1:S
    SING=[SING,double(RSIN{j}+RSID{j})];
end
AA=SING;


%% compute contributions to other sides

r1 = @(t) sqrt(1+(1+t).^2/4);
r2 = @(t) sqrt(1+t.^2);

theta1 = @(t) atan((1+t)/2);
theta2 = @(t) acos(-t./(sqrt(1+t.^2)));%atan(-1./t) + (t>=0)*pi;

x = chebfun('x');
[RD1,~,~,~] = AB(z(1),z(2),k,50);
[~,RN2,~,~] = AB(z(2),z(3),k,50);

for j=1:S
    j
    Q = @(t) r1(t).^(2/3*(2*j-1)).*sin((2/3)*(2*j-1)*theta1(t));  
    p=Q(x);
    a = transpose(chebpoly(p));
    a = a(end:-1:1);
    b1=cheb2leg(a);
    length(b1);    
    
    AA(:,j)=AA(:,j)+RD1(:,1:length(b1))*b1(1:length(b1));
    
    Q = @(t) r2(t).^(2/3*(2*j-1)-1)*(2/3)*(2*j-1).*(sin(theta2(t)).*sin((2/3)*(2*j-1)*theta2(t))+cos(theta2(t)).*cos((2/3)*(2*j-1)*theta2(t)));
    p=Q(x);
    a = transpose(chebpoly(p));
    a = a(end:-1:1);
    b1=cheb2leg(a);
    length(b1);

    AA(:,j)=AA(:,j)+RN2(:,1:length(b1))*b1(1:length(b1));
    
end
    
    
    


%%

AA=[AA;conj(AA)];
BB=[BB;conj(BB)];
A1=AA;
B1=BB;

qd2=zeros(N,1);
qn1=qd2;
qd2(1)=1;
rB=[qn1;qd2];

    
%pre-condition

sc = 1./sum(abs([AA,BB]),2);
AA  = bsxfun(@times,sc,AA);
RH = BB*rB;
RH  = bsxfun(@times,sc,RH);

sol = real(AA\RH);
cjl=sol(1:end);
end
