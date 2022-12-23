%|Numerical DFA Code|University of Illinois at Urbana-Champaign| 
%|ME 370|Created by Brian C. McGuigan

%Section1:Input
cmarray=[0.5;0.5;0.5;0.5;;0.5;0.5;0.5;0.5;0.5;0.5];                  %Center of mass array. Single column array that holds a value for each link (0-1) which represents the fractional link distance from the first node entry
                           % that the link's center of mass is. Supply an
                           % entry for the ground link as well even though it has
                           % no moving center of mass
gconst=9.81; %gravitational constant in m/s^2 (positive)
wdth=1*0.0254; %Width of each link (meters). (0.0254m=1inch)
tks=1*0.003175; %thickness of each link (meters). (1/8inch=0.003175m)
den=1180; %Density for each link (kg/m^3). (acrylic: 1180 kg/m^3)
Mass=0.624;
coeff_fric=0.0;
%Section2: Calculates the mass of each link using the previously defined link lengths
%and the above constant width, thickness, and density
M=zeros(length(cmat),1); %Allocate single column array for mass of each link 
for i=1:length(cmat)
    M(i)=den*links(i)*wdth*tks; %calculates mass of link i
end

%Section3: Calculates the moment of inertia for each link assuming the above
%rectangular link dimension and relative center of mass position from
%cmarray
I= zeros(length(cmat),1); %Moment of inertias for each link. Single column matrix where row i is the moment of inertia of link i
for i=1:length(cmat)
   
   I(i)=0.5*M(i)*(wdth^2+links(i)^2); %moment of inertia of rectangular link through center
   I(i)=I(i)+M(i)*((0.5-cmarray(i))*links(i))^2; %modify according to center of mass location using parallel axis theorem. If cmarray(i)=0.5, no change is made because the center of mass is through the center
end

%Section4: Calculates the angular velocity of each link
omega=zeros(length(inodes),tnum-1); %omega for each link in radians
for k=1:tnum-1
    for i=1:length(cmat)
        nd1ind=cmat(i,1); %node 1 index of link i
        nd2ind=cmat(i,2); %node 1 index of link i
        v1=[real(velocity(nd1ind,k)) imag(velocity(nd1ind,k)) 0]; %velocity of node 1 in x,y,z array form
        v2=[real(velocity(nd2ind,k)) imag(velocity(nd2ind,k)) 0]; %velocity of node 2 in x,y,z array form
        r1=[real(xnode(nd1ind,k)) imag(xnode(nd1ind,k)) 0]; %position of node 1 in x,y,z array form
        r2=[real(xnode(nd2ind,k)) imag(xnode(nd2ind,k)) 0]; %position of node 2 in x,y,z array form
        w=cross(v1-v2,r2-r1)/norm(r2-r1)^2; %calculation of omega in array form
        omega(i,k)=w(3); %gets z component which is axis of rotation for our 2d constrained
    end
end

%Section5: Calculates the angular acceleration of each link
alpha=zeros(length(inodes),tnum-2); %alpha for each link in radians
for k=1:tnum-2
    for i=1:length(cmat)
        nd1ind=cmat(i,1); %node 1 index of link i
        nd2ind=cmat(i,2); %node 2 index of link i
        a1=[real(acceleration(nd1ind,k)) imag(acceleration(nd1ind,k)) 0]; %acceleration of node 1 in x,y,z array form
        a2=[real(acceleration(nd2ind,k)) imag(acceleration(nd2ind,k)) 0]; %acceleration of node 2 in x,y,z array form
        r1=[real(xnode(nd1ind,k)) imag(xnode(nd1ind,k)) 0]; %position of node 1 in x,y,z array form
        r2=[real(xnode(nd2ind,k)) imag(xnode(nd2ind,k)) 0]; %position of node 2 in x,y,z array form
        a=cross(a1-a2+norm(omega(i,k))^2*(r1-r2),r2-r1)/norm(r2-r1)^2; %calculation of alpha in array form
        alpha(i,k)=a(3); %gets z component which is axis of rotation for our 2d constrained
    end
end

%Section6: Calculates the linear acceleration of each link's of
%mass/gravity by going from the acceleration of the first node of each link
%to the acceleration of the center of gravity using calculated alpha of
%each link
acg=zeros(length(inodes),tnum-2); %linear acceleration for each link
for k=1:tnum-2 %loops over all timesteps of simulation
    for i=1:length(cmat) %loops over all links
        
        nd1ind=cmat(i,1); %node 1 index of link i
        nd2ind=cmat(i,2); %node 2 index of link i
        r2=xnode(nd2ind,k); %Position of node 2 for link
        r1=xnode(nd1ind,k); %Position of node 1 for link
        
        rcg=cmarray(i)*(r2-r1)+r1; %Position of center of gravity for link as complex number.
        a1=[real(acceleration(nd1ind,k)) imag(acceleration(nd1ind,k)) 0]; %linear acceleration at node 1
        r1cg=[real(rcg-r1) imag(rcg-r1) 0]; %Distance from r1 to rcg
        
        acgcomp=a1+cross([0 0 alpha(i,k)],r1cg)-norm(omega(i,k))^2*r1cg; %Acg in x,y,z array form
        acg(i,k)=acgcomp(1)+1i*acgcomp(2); %puts x,y,z array form to complex notation
        
    end
end

%Section7: Allocates B
nlinkeqs=3*(length(cmat)-1); %Number of link equations. Sum of Fx,Fy,Torque for all links except ground.
nnodeeqs=2*length(inodes); %Number of nodal equations. Sum of Fx, Fy at each node. 
B=zeros(nlinkeqs+nnodeeqs-2,tnum-2);%Allocate solution vector


%Makes an ordered array of all link indices not including the ground link.
%lnksnotgnd is used when filling in [A] matrix
lnksnotgnd=1:length(cmat); %makes ordered array of all links in system
len=length(gndlink);
lnksnotgnd=lnksnotgnd(lnksnotgnd~=gndlink(1)); %removes gndlink index

for j=1:tnum-2 %Loops through each timestep
    
    %Allocates and rezeros A,C
    A=zeros(nlinkeqs+nnodeeqs-2,nlinkeqs+nnodeeqs-2);
    C=zeros(nlinkeqs+nnodeeqs-2,1);

    %Section7: Fills in link components (Sum Fx,Fy,M) in A,C matrix/vector
    for i=1:length(cmat)-1 %loops over rows of A matrix by link index
        idx=lnksnotgnd(i);
        
        %Fills in Sum of Fx matrix/equation lines for link i
            %Insert fx coefficients on A
            
            A(3*(i-1)+1,(idx-1)*4+1)=1;
            A(3*(i-1)+1,(idx-1)*4+3)=1;
            
            %Insert in mass and x-acceleration on C
            
            C(3*(i-1)+1)=M(idx)*real(acg(idx,j));
          
        %Fills in Sum of Fy matrix/equation lines for link i
            %Insert fy coefficients on A
            A(3*(i-1)+2,(idx-1)*4+2)=1;
            A(3*(i-1)+2,(idx-1)*4+4)=1;
           
            %Insert in mass and y-acceleration on C, should include graviational force
            
            C(3*(i-1)+2)=M(idx)*(imag(acg(idx,j))+gconst);
          
        %Fills in Sum of moment matrix/equation lines for link i
        
        r1=xnode(cmat(idx,1),j);    %Insert position of node 2 for link i
        r2=xnode(cmat(idx,2),j);    %Insert position of node 1 for link i
        rcg=cmarray(lnksnotgnd(i))*(r2-r1)+r1; %Position of center of gravity of link i as complex number (excluding ground)
        r1cg=r1-rcg;
        r2cg=r2-rcg;
        
            %Insert first moment equation coefficient for node 1
            A(3*i,(idx-1)*4+1)=-imag(r1cg);
                       
            %Insert second moment equation coefficient for node 1
            A(3*i,(idx-1)*4+2)=real(r1cg);
            
            %Insert first moment equation coefficient for node 2
            A(3*i,(idx-1)*4+3)=-imag(r2cg);
            
            %Insert second moment equation coefficient for node 2
            A(3*i,(idx-1)*4+4)=real(r2cg);
            
            %Insert moment of inertia and angular accel component in C
            C(3*i)=I(idx)*alpha(idx,j);
        
        %Fills in moment coefficient in matrix for the motor node
        if lnksnotgnd(i)==crlink
            %Inser moment coefficient for motor node on crank link
            A(3*i,nlinkeqs+nnodeeqs-2)=1;
        end
        
    end
    
    %Section8: Fills in node components to stiffness matrix
    for i=1:length(inodes) %loops over node indices to fill in matrix components
        
        [linkind,nodecol]=find(cmat==i); %finds link indices that share the ith node and the column the node is in
        for k=1:length(linkind) %loops over links that share node i
            
            %Insert x component of node action-reaction equations
            A(nlinkeqs+2*(i-1)+1,4*(linkind(k)-1)+2*(nodecol(k)-1)+1)=1;
            
            %Insert y component of node action-reaction equations
            A(nlinkeqs+2*(i-1)+2,4*(linkind(k)-1)+2*(nodecol(k)-1)+2)=1;
        end
        factor=0;
        if(crankangles(j)>180)
            factor=1;
        end
        if(i==1)
            C(nlinkeqs+2*(i-1)+1)=sign(crangvel)*Mass*gconst*coeff_fric*factor;
            C(nlinkeqs+2*(i-1)+2)=Mass*gconst*factor*2.2;
        else
            C(nlinkeqs+2*(i-1)+1)=0;
            C(nlinkeqs+2*(i-1)+2)=0;
        end    
    end
    
    %Section9: Solves for force vector unknowns.
    B(:,j)=pinv(A)*C;
end

Torque=B(nlinkeqs+nnodeeqs-2,:);

t_plot=figure();%=figure(k_x)

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(crankangles,Torque);

hxlabel=xlabel('Angle(deg)');
hylabel=ylabel('Torque(Nm)');
htitle=title('Torque Plot');
formatplot(t_plot,hax,hxlabel,hylabel,htitle);
legend('Torque');

%Links 6 9 3 1

% idx=1;
% 
% r_plot=figure(2);%=figure(k_x)
% 
% hax=axis(); %initiate axis
% hold on
% 
% R1=abs(B((idx-1)*4+1,:)+1i*B((idx-1)*4+2,:));
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(crankangles,R1);
% 
% hxlabel=xlabel('Angle(deg)');
% hylabel=ylabel('Force(N)');
% htitle=title('Reaction at Node 1');
% formatplot(r_plot,hax,hxlabel,hylabel,htitle);
% legend('Force');
% 
% idx=1;
% 
% q_plot=figure(3);%=figure(k_x)
% 
% hax=axis(); %initiate axis
% hold on
% 
% R2=abs(B((idx-1)*4+3,:)+1i*B((idx-1)*4+4,:));
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(crankangles,R2);
% 
% hxlabel=xlabel('Angle(deg)');
% hylabel=ylabel('Force(N)');
% htitle=title('Reaction at Node 2');
% formatplot(q_plot,hax,hxlabel,hylabel,htitle);
% legend('Force');
% 
% idx=3;
% 
% p_plot=figure(4);%=figure(k_x)
% 
% hax=axis(); %initiate axis
% hold on
% 
% R3=abs(B((idx-1)*4+1,:)+1i*B((idx-1)*4+2,:));
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(crankangles,R3);
% 
% hxlabel=xlabel('Angle(deg)');
% hylabel=ylabel('Force(N)');
% htitle=title('Reaction at Node 3');
% formatplot(p_plot,hax,hxlabel,hylabel,htitle);
% legend('Force');
% 
% idx=3;
% 
% o_plot=figure(5);%=figure(k_x)
% 
% hax=axis(); %initiate axis
% hold on
% 
% R4=abs(B((idx-1)*4+3,:)+1i*B((idx-1)*4+4,:));
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(crankangles,R4);
% 
% hxlabel=xlabel('Angle(deg)');
% hylabel=ylabel('Force(N)');
% htitle=title('Reaction at Node 4');
% formatplot(o_plot,hax,hxlabel,hylabel,htitle);
% legend('Force');
% 
% idx=6;
% 
% m_plot=figure(6);%=figure(k_x)
% 
% hax=axis(); %initiate axis
% hold on
% 
% R5=abs(B((idx-1)*4+1,:)+1i*B((idx-1)*4+2,:));
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(crankangles,R5);
% 
% hxlabel=xlabel('Angle(deg)');
% hylabel=ylabel('Force(N)');
% htitle=title('Reaction at Node 5');
% formatplot(m_plot,hax,hxlabel,hylabel,htitle);
% legend('Force');
% 
% idx=6;
% 
% k_plot=figure(7);%=figure(k_x)
% 
% hax=axis(); %initiate axis
% hold on
% 
% R6=abs(B((idx-1)*4+3,:)+1i*B((idx-1)*4+4,:));
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(crankangles,R6);
% 
% hxlabel=xlabel('Angle(deg)');
% hylabel=ylabel('Force(N)');
% htitle=title('Reaction at Node 6');
% formatplot(k_plot,hax,hxlabel,hylabel,htitle);
% legend('Force');
% 
% idx=9;
% 
% b_plot=figure(8);%=figure(k_x)
% 
% hax=axis(); %initiate axis
% hold on
% 
% R7=abs(B((idx-1)*4+1,:)+1i*B((idx-1)*4+2,:));
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(crankangles,R7);
% 
% hxlabel=xlabel('Angle(deg)');
% hylabel=ylabel('Force(N)');
% htitle=title('Reaction at Node 7');
% formatplot(b_plot,hax,hxlabel,hylabel,htitle);
% legend('Force');
% 
% idx=9;
% 
% t_plot=figure(9);%=figure(k_x)
% 
% hax=axis(); %initiate axis
% hold on
% 
% R8=abs(B((idx-1)*4+3,:)+1i*B((idx-1)*4+4,:));
% 
% hax = gca; %initialize current axis
% hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
% plot(crankangles,R8);
% 
% hxlabel=xlabel('Angle(deg)');
% hylabel=ylabel('Force(N)');
% htitle=title('Reaction at Node 8');
% formatplot(t_plot,hax,hxlabel,hylabel,htitle);
% legend('Force');