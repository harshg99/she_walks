close all
clear
clc
% Created by Harsh Goel to compile the results of DFA
%Variable for tracking different graphs
k_x=1;
%stroing maximum torques at different rotation speeds
MaxT=zeros(9,1);
for j=6:4:30
    %angular velocity of crank from rpm to deg/s
    crangvel=-j*6;
    
    %running analysis
    PVA;
    k_x=k_x+1;
    DFA_main;
    
    %stroing maximum torque
    MaxT(k_x)=max(Torque);
end

%plotting maximum torque
hfig=figure(10);
hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot([6,10,14,18,22,26,30],MaxT(2:8));
    
hxlabel=xlabel('Angular Velocity (rpm)');
hylabel=ylabel('Max Torque(Nm)');
htitle=title('Max Torque Plot');
formatplot(t_plot,hax,hxlabel,hylabel,htitle);
legend('Max Torque');

crangvel=10*6;
% PVA;
% k_x=12;
% DFA_main;
% 
% %Energy change verification
% Pin=zeros(tnum-2);
% Pout=zeros(tnum-2);
% 
% vel=zeros(length(inodes),tnum-1); %omega for each link in radians
% for k=1:tnum-1
%     for i=1:length(cmat)
%         nd1ind=cmat(i,1); %node 1 index of link i
%         nd2ind=cmat(i,2); %node 1 index of link i
%         v1=[real(velocity(nd1ind,k)) imag(velocity(nd1ind,k)) 0]; %velocity of node 1 in x,y,z array form
%         v2=[real(velocity(nd2ind,k)) imag(velocity(nd2ind,k)) 0]; %velocity of node 2 in x,y,z array form
%         vel(i,k)=velocity(nd1ind,k)/2+velocity(nd2ind,k)/2;
%         
%     end
% end
% for j=1:tnum-2
%     %power input in W
%     Pin(j)=Torque(j)*10/60*2*pi;
%     %net power output in W
%     for i=1:length(cmat)-1
%         Pout(j)=Pout(j)+M(idx)*(gconst*imag(vel(i,j))+real(acg(i,j))*real(vel(i,j))+imag(acg(i,j))*imag(vel(i,j)))+I(idx)*omega(i,j)*alpha(i,j);        
%     end
% end
% 
% hfig=figure(13);
% hax=axis(); %initiate axis
% hold on
% 
% 
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(dt:dt:(tnum-2)*dt,Pin);
% plot(dt:dt:(tnum-2)*dt,Pout);
%     
% hxlabel=xlabel('Time(s)');
% hylabel=ylabel('Power/Rate of Change of Energy(W)');
% htitle=title('Power Plot');
% formatplot(t_plot,hax,hxlabel,hylabel,htitle);
% legend('Power Input','Power Output');


    