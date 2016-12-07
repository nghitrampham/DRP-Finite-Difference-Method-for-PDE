%heat2D_driver
close all;
Tf = (.1:0.1:4);
%Tf = (1:.2:3);
figure(1);

for i = 1:length(Tf)
    %[x ,y, uApproxTf] = heat2D(@(x,y,t) sin(pi*x*t)*cos(pi*y*t), ... 
     %       @(x,y) cos(pi*y + (pi/2))*cos(pi*x + (pi/2)), @(x,y,t) 0, 1, Tf(i), 15);
      [x ,y, uApproxTf] = heat2D(@(x,y,t) 0  , ... 
             @(x,y) (-1/exp(1)+exp(-(2*x-1).^2)).*(-1/exp(1)+exp(-(2*y-1).^2)), @(x,y,t) 0, 1, Tf(i), 15);   
    %figure(i);
    surf(x,y,uApproxTf);
    
    xlim([0 1]);
    ylim([0 1]);
    zlim([-.04 .04]);
    drawnow;
    pause(.5);
end