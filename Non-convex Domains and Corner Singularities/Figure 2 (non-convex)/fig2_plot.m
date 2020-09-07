close all

Col.R=50;
Col.K=180;
Col.collocation_type='halton';

tvec=0:0.01:.95;

[ max_err1 ] = max_error1( Col,tvec ,14);
[ max_err2 ] = max_error2( Col,tvec ,14);

Col.collocation_type='ray';
Col.R2=30;
Col.R1=Col.R2/(Col.K/4);

[ max_err3 ] = max_error1( Col,tvec,14 );
Col.K=180;
Col.collocation_type='ray';

Col.R1=Col.R2/(Col.K/3);
[ max_err4 ] = max_error2( Col,tvec,14 );
%%
figure
FS = 15;
semilogy(tvec,max_err1,'linestyle','--','linewidth',2,'Color','b');% [0.7 .7 .7])
hold on
semilogy(tvec,max_err2,'linewidth',2,'Color','b')% [0.7 .7 .7])
semilogy(tvec,max_err3,'--k')
semilogy(tvec,max_err4,'k')
plot([1/3,1/3],[10^(-16),100],'r','linewidth',1)

ax = gca;
set(ax,'YTick',10.^(-16:2));
% axis([0,1,10^-10,10])

xlabel('$x_4$','interpreter','latex','FontSize',FS);
ylabel('Error (log-scale)','interpreter','latex','FontSize',FS);

legend({'Halton without internal edge','Halton with internal edge','Ray without internal edge','Ray with internal edge','Convexity Transition'},'interpreter','latex','FontSize',FS,'location','Northwest')
axis([0,1,10^(-13),10])