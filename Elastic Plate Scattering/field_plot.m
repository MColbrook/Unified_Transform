%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SET PARAMETERS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Physical parameters
close all
clear
K0=10; % freq of incidnet field
theta=pi/4; % angle of incident field
pos=[0,1]; % endpoints of plates/gaps
ty=[4]; % type of plates/gaps, 0=gap,1=stiff,2=CC,3=FF,4=CF,5=FC

% physical parameters for elastic plates
epsilon=0.0021;
Omega=[0.25];
k=K0^4./Omega.^4;
omega=K0;
B=Omega.^6/(epsilon*K0^3);
% the following are vectors to allow plates with different properties
C1=0*ty+B./(omega.^2); % first constant defining properties of plate, input as a vector if this is different for each plate
C2=-k.*C1; % second constant defining properties of plate, input as a vector if this is different for each plate
C3=0*C1+omega.^2; % this is used to compute the plate displacement

%% Numerical parameters

N=zeros(1,length(pos)+1)+100; % number of basis function in each segment (including infinite parts) from left to right
X=-2.5:0.01:2.5; % these form a grid at points where we want to plot the field
Y=-2:0.01:2; % these form a grid at points where we want to plot the field
L1=length(X);
L2=length(Y);
[X0,Y0] = meshgrid(X,Y);
Z=X0(:)+1i*Y0(:); % these form a grid at points where we want to plot the field

x=cell(length(pos)-1,1); % these are points on each plate where we compute the jump in the scattered field
x{1}=-1:0.01:1; % scaled to interval [-1,1] for ease
for j=2:length(pos)-1
    x{j}=x{1};
end
y=x; % these are points on each plate where we compute the plate displacement

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COMPUTATION %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Q_ref,~,~,q_ref,eta_ref,sol_ref] = UT_multiplate_collinear(K0,theta,N,pos,ty,C1,C2,C3,Z,[],x,y);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOT RESULTS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I=[];
for j=1:length(pos)-1
    if ty(j)>0
        I=union(find((abs(imag(Z))<0.015)&(real(Z)<(.015+pos(j+1)))&(real(Z)>(pos(j)-0.01))),I(:));
    end
end
q0=exp(-1i*K0*(cos(theta)*real(Z)+sin(theta)*imag(Z))); % incident field
q=Q_ref+q0; % obviously you can change this to just the scattered field
q(I)=NaN; % efficient way to plot the plate
Q_ref(I)=NaN;
q=reshape(q,L2,L1);
Q_ref=reshape(Q_ref,L2,L1);
map=zeros(100,3);
for j=1:100
    map(j,1)=(j-1)/99;
    map(j,3)=1-(j-1)/99;
    map(j,:)=[1,1,1]*(1-abs(j-50.5)/50)+map(j,:)*(abs(j-50.5)/50);
end

dp = [-cos(theta),-sin(theta)]/2;
pp1=[0.5*(pos(end)+pos(end-1)),0]-2*dp;

figure
imagesc(X,Y,real(q),'AlphaData',~isnan(q));
set(gca,'YDir','normal')
colormap(map)
hold on
title('Real Part of Total Field')
colorbar
set(gca,'color',0*[1 1 1]);
quiver(pp1(1),pp1(2),dp(1),dp(2),0,'k','linewidth',1,'MaxHeadSize',1)
axis equal tight

figure
imagesc(X,Y,abs(q),'AlphaData',~isnan(q));
set(gca,'YDir','normal')
colormap(map)
hold on
title('Absolute Value of Total Field')
colorbar
set(gca,'color',0*[1 1 1]);
quiver(pp1(1),pp1(2),dp(1),dp(2),0,'k','linewidth',1,'MaxHeadSize',1)
axis equal tight

figure
imagesc(X,Y,real(Q_ref),'AlphaData',~isnan(q));
set(gca,'YDir','normal')
colormap(map)
hold on
title('Real Part of Scattered Field')
colorbar
set(gca,'color',0*[1 1 1]);
quiver(pp1(1),pp1(2),dp(1),dp(2),0,'k','linewidth',1,'MaxHeadSize',1)
axis equal tight

figure
imagesc(X,Y,abs(Q_ref),'AlphaData',~isnan(q));
set(gca,'YDir','normal')
colormap(map)
hold on
title('Absolute Value of Scattered Field')
colorbar
set(gca,'color',0*[1 1 1]);
quiver(pp1(1),pp1(2),dp(1),dp(2),0,'k','linewidth',1,'MaxHeadSize',1)
axis equal tight
%% Uncomment to plot plate deformations
 
% figure
% for j=1:(length(pos)-1)
%     if ty(j)>0
%         y1=real(eta_ref{j});
%         x1=(y{j}+1)/2*(pos(j+1)-pos(j))+pos(j);
%         plot(x1,y1,'k')
%         hold on
%     end
% end
% axis([min(X),max(X),-0.1,0.1])











