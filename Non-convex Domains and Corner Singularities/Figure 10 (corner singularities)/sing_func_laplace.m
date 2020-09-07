clc
clear
close all


N2vec=5:29;
L=length(N2vec);

for j=1:L
    j
    N=N2vec(j);
    col.type='ray';
    col.M=2*N;
    col.R1=1/10;
    col.R2=2*N;
    col.K=4*N;
    col.R=10;
    
    [ cjl,D] = Lap_sing( 1,col,N2vec(j));
    cjl;
    
    Out1{j}=cjl;
    C1(j)=cond(D);
    
    %% current error
    max(abs(cjl(1:5)'-[1.127980401059388,0.169933866502253,-0.023040973993480,0.003471196658216,0.000915157099087]))

end



%%
figure
clear E1 E2 E3 E4 E5
for j=1:L
    E1(j)=Out1{j}(1);
    E2(j)=Out1{j}(2);
    E3(j)=Out1{j}(3);
    E4(j)=Out1{j}(4);
    E5(j)=Out1{j}(5);
end

V=[1.127980401059388;  % values computed using larger N and testing against other collocation points
     0.169933866502253;
     -0.023040973993480;
     0.003471196658216;
     0.000915157099087;
     -0.000112803834456;
     0.000087716524487;
     0.000027760313669;
     -0.000004416157802;
     0.000002753945681;
     0.000000921961930;
     -0.000000155445876;
     0.000000108840774;
     0.000000037969831;
     -0.000000006661925;
     0.000000004710615;
     0.000000001682649;
     -0.000000000301797;
     0.000000000218392;
     0.000000000079333;
     -0.000000000014456;
     0.000000000010590;
     0.000000000003894;
     -0.000000000000718;
     0.000000000000531;
     0.000000000000197;
     -0.000000000000037;
     0.000000000000027;
     0.000000000000010];
 
for j=1:L
    V1=Out1{j};
    V1=V1(:);
    E(j)=max(abs(V(1:length(V1))-V1));
end

%% Error in interior points
load('points.mat')

Z_sol=0*r;
for mu=1:length(V)
    Z_sol=Z_sol+V(mu)*r.^(2/3*(2*mu-1)).*sin(2/3*(2*mu-1)*theta);
end

for j=1:L
    V1=Out1{j};
    Z=0*r;
    for mu=1:length(V1)
        Z=Z+V1(mu)*r.^(2/3*(2*mu-1)).*sin(2/3*(2*mu-1)*theta);
    end
    E_inside(j)=max(abs(Z-Z_sol));
end

%%

semilogy(N2vec(1:L),abs(E1-1.127980401059388),'-ok');
hold on
semilogy(N2vec(1:L),abs(E2-0.169933866502253),'-sk');
semilogy(N2vec(1:L),abs(E3+0.023040973993480),'-^k');
semilogy(N2vec(1:L),abs(E4-0.003471196658216),'-*k');
semilogy(N2vec(1:L),abs(E5-0.000915157099087),'-xk');


semilogy(N2vec(1:L),E,'--k');
semilogy(N2vec,E_inside,'b','linewidth',2)
axis([4,30,10^(-16),10^(-3)])

xlabel('Number of Basis Functions $N$','interpreter','latex','fontsize',15)
ylabel('Error','interpreter','latex','fontsize',15)
legend({'$\alpha_1$','$\alpha_2$','$\alpha_3$','$\alpha_4$','$\alpha_5$','$l^\infty$','Interior'},'interpreter','latex','fontsize',12,'location','NorthEast')


