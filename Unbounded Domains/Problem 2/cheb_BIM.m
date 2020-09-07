function [q]=cheb_BIM(K0,theta,N)
%%% compute approximation via boundary integral method

K=zeros(N);

fun1 = @(t,m,n) (sqrt(1-1./t.^2)-1).*besselj(n,K0*t/2).*besselj(m,K0*t/2)./t;
fun2 = @(t,m,n) (sqrt(1-1./t.^2)-1+1./(2*t.^2)).*besselj(n,0.5*K0*t).*besselj(m,t*K0/2)./t;
fun3 = @(t,m,n) (sqrt(1-1./t.^2)-1+1./(2*t.^2)+1./(8*t.^4)).*besselj(n,0.5*K0*t).*besselj(m,0.5*K0*t)./t;

for aa=1:N
    for bb=1:N
        if mod(aa+bb,2)==1
        elseif bb>aa-1
            if aa+bb >5
                K(aa,bb)=quadgk(@(t) fun3(t,aa,bb),0,Inf,'AbsTol',eps,'MaxIntervalCount',10^4)+(aa==bb)/(2*aa);
                    if aa==bb
                        K(aa,bb)=K(aa,bb)-K0^2/32/(aa*(aa-1)*(aa+1))-K0^4/(16*64*aa*(aa+1)*(aa+2))*3/(2*(aa-1)*(aa-2));
                    elseif aa==(bb-2)
                        K(aa,bb)=K(aa,bb)-K0^2/32/(2*aa*(aa+2)*(aa+1))-K0^4/(16*64*aa*(aa+1)*(aa+2))/((aa+3)*(aa-1));
                    elseif aa==(bb-4)
                        K(aa,bb)=K(aa,bb)-K0^4/(16*64*aa*(aa+1)*(aa+2))/(4*(aa+3)*(aa+4));
                    end                 
            elseif aa+bb>3
                K(aa,bb)=quadgk(@(t) fun2(t,aa,bb),0,Inf,'AbsTol',eps,'MaxIntervalCount',10^4)+(aa==bb)/(2*aa);
                if aa==bb
                    K(aa,bb)=K(aa,bb)-K0^2/32/(aa*(aa-1)*(aa+1));
                elseif aa==(bb-2)
                    K(aa,bb)=K(aa,bb)-K0^2/32/(2*aa*(aa+2)*(aa+1));
                end
            else
                K(aa,bb)=quadgk(@(t) fun1(t,aa,bb),0,Inf,'AbsTol',eps,'MaxIntervalCount',10^4)+(aa==bb)/(2*aa);
            end
        end
    end
    K(aa,aa)=K(aa,aa)/2;
end

K=(K+transpose(K));
K=conj(K); % due to def of branch cut
b=zeros(N,1);
for j=1:N
    b(j)=-1i*besselj(j,K0/2*cos(theta))*tan(theta);
end
sol=transpose((K\b)*2)./(1:N);


t=-1:0.01:1;
q=0*t;
for j=1:N
    q=q+sol(j)*Chebval(j,t);
end
end


function [X]=Chebval(n,z)
    if mod(n,2)==0
        X=1i*sin(n*asin(z));
    else
        X=cos(n*asin(z));
    end
end







                
            