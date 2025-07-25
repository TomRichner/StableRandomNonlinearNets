function plot_hhrk4(tt,yy,stim)

clf
plot(tt,stim(1:2:(2*length(tt))));
title('Stimulus')

figure
subplot(2,1,1)
plot(tt,yy(1,:))
title('HH mV vs msec')

subplot(2,1,2)
plot(tt,yy(2,:),'b',tt,yy(3,:),'g',tt,yy(4,:),'r')
legend('n','m','h');

