%|Numerical PVA Code|University of Illinois at Urbana-Champaign| 
%|ME 370|Created by Brian C. McGuigan
close all
clear
clc

% Section1:Input

inodes=0.0254*[-0.0015 + 0.6047i;
        1.0879+2.1918i;
        1.8762+2.9551i;
        2.769+3.1092i;
        4.1395+2.6641i;
        3.6033+2.5415i;
        2.769+1.9705i;
        2.1501+2.2089i];     % Node/Joint Placement
start=inodes(1);
inodes=inodes-start;
scale=4/3.18;
inodes=inodes*scale;
cmat=[1,2;
      2,3;
      3,4;
      4,5;
      4,7;
      5,6;
      6,8;
      5,7;
      7,8;
      8,2];       % Connectivity Matrix
gndlink=[4,5,8];	% Link Index of Ground Link (row of cmat denoting nodes to freeze during simulation)
crlink=6;     % Link with prescribed position over time (input crank)
motornode=5;	% Node the motor is at
slider={};  % Slider constraints
fxdangle=[1,2;7,10];% Fix link/link angle during rotation. (Must share a node)
rotfix=[];  % Nodes that are fixed so that they cant rotate
slfric=[];  % Friction at slider constraints

link_length=zeros(10,1);
for j=1:10
    link_length(j)=200*norm(inodes(cmat(j,1))-inodes(cmat(j,2)));
end

tperiod=2.5;	% total time to run simulation (seconds)
dt=0.1;         % time step (seconds)
tnum=tperiod/dt; %Total number of timesteps
time=0:dt:tperiod; %array holding the amount of time that has passed
crangvel=47.5/60*360;   % Angular velocity of input crank (deg/sec) -Assumed constant.
             %+(-) angular velocity denotes CCW (CW) rotation

% Section2: Useful information about the sytem
%fxnodes=unique(setdiff(cmat([crlink,gndlink],:),motornode)); %Nodes with prescribed/fixed placement (the ones that make up the ground and crank link)
fxnodes=unique(cmat([gndlink],:));
 if ~isempty(slider)
     for j=1:length(slider)
        fxnodes=setdiff(fxnodes,slider{j}(1));
     end
end
 
nfxnodes=setdiff(1:length(inodes),fxnodes); %Nodes without prescribed/fixed placement (the ones that can vary after we move crank link a bit)
%crpvtnd= intersect(cmat(gndlink,:),cmat(crlink,:));%Node index that input crank link rotates about
crpvtnd=motornode; %Assigns previous variable entry
crtipnd= setdiff(cmat(crlink,:),crpvtnd);%Node index denoting the tip of the input crank link
ndof=length(nfxnodes)*2; %Nu(3*i,end)=1; %Torque coeefficient on crank link
linkind=setdiff(1:length(cmat),gndlink);%links that may have varying lengths after a timestep

% Slider information used for DFA
sliderlinks=[];
sliderendnodes=[];
if ~isempty(slider)
     for j=1:length(slider)
         linkind=setdiff(linkind,slider{j}(2));
        if length(slider{j})==3 
            linkind=setdiff(linkind,slider{j}(3));
            sliderlinks=[sliderlinks;j, slider{j}(2), slider{j}(3);];%total list of all slider links [sliderindex, link1ind linke1ind]
            cmnode=intersect(cmat(slider{j}(2),:),cmat(slider{j}(3),:));
            endnodes=setdiff([cmat(slider{1}(2),:) cmat(slider{1}(3),:)],cmnode);
            sliderendnodes=[sliderendnodes;endnodes(1),endnodes(2);];
        end
     end
end

% Section3: Allocation/Initialization
xnode=zeros(length(inodes),tnum); %allocate node positions over all time steps 
xnode(:,1)=inodes; %Initialize current timestep (put in initial node position data)

% Section4: Put in initial fixed ground link nodes positions at each time step in corresponding position array since they are fixed
for i=2:tnum
    for j=1:length(gndlink)
        xnode(cmat(gndlink(j),:),i)=inodes(cmat(gndlink(j),:));
    end
end


% Section6: Computes link lengths from initial node positions
links=zeros(length(cmat),1); %Allocate array for link lengths
for i=1:length(cmat) %loop over connectivity matrix
    % Fill in the remaining code here
    l1=cmat(i,:);
    links(i)=abs(inodes(l1(1))-inodes(l1(2)));
end

%Compute initial link angle for fixed angle constraint
angs=zeros(size(fxdangle,1),1);
angshnode=zeros(size(fxdangle,1),1);
rotangs=cell(size(rotfix,1),1);

%Calculates the initial angles for the fixed angle constraints
if ~isempty(angs)
    for i=1:size(fxdangle,1)
            
        ndind1=intersect(cmat(fxdangle(i,1),:),cmat(fxdangle(i,2),:));
        ndind2=setdiff(cmat(fxdangle(i,1),:),ndind1);
        vec1=inodes(ndind2)-inodes(ndind1);
      
        ndind2=setdiff(cmat(fxdangle(i,2),:),ndind1);
        vec2=inodes(ndind2)-inodes(ndind1);      
        angs(i)=(real(vec1)*real(vec2)+imag(vec1)*imag(vec2))/(norm(vec1)*norm(vec2));
    end
end

%Calculates the initial link vectors from each rotation fixed constraint location
if ~isempty(rotangs)
    for i=1:length(rotfix)        
        [linkrow,nodecol]=find(cmat==rotfix(i));
        rotangs{i}=zeros(length(linkrow),2);
        for k=1:length(linkrow)
            node2in=cmat(linkrow(k),setdiff(1:2,nodecol(k)));
            node1in=rotfix(i);
            vec=inodes(node2in)-inodes(node1in);
            
            rotangs{i}(k,1)=real(vec);
            rotangs{i}(k,2)=imag(vec);
        end
        
    end
end
%Section7: Loop over each time step and find adjustments that preserve link lengths
str=sprintf('Running Simulation:'); %Make string
disp(str) %Display String
linksnew=zeros(length(linkind),1); %Allocate array for new link lengths (link lengths that could change in objective function after dx is solved
for i=1:tnum-1
    
    xnode(nfxnodes,i+1)=xnode(nfxnodes,i);%copy previous nfxnode quantities to current timestep
    options = optimset('Display', 'off'); %Turns off display options
    %options.OptimalityTolerance=1e-12;
    %options.FunctionTolerance=1e-12;
    %options.StepTolerance=1e-12;
    %options.Algorithm = 'levenberg-marquardt';
    func=@(dx)linklen(dx,dt,xnode(:,i+1),linkind,cmat,links(:),nfxnodes,slider,fxdangle,angs,crangvel,inodes,motornode,crtipnd,rotfix,rotangs); %Sets function handle to include additional arguments
    dxnfx=lsqnonlin(func,1e-6*ones(1,ndof),[],[],options); %Solve for nfxnode perturbation distance that fixes link length incompatibilities
    xnode(nfxnodes,i+1)=xnode(nfxnodes,i+1)+dxnfx(1:(ndof/2))'; %add real part correction
    xnode(nfxnodes,i+1)=xnode(nfxnodes,i+1)+dxnfx((ndof/2+1):end)'*1i; % add imag part correction
    %motornodevel=xnode(motornode,i+1)-xnode(motornode,i);
    
    %Loops over each link that can have varied length to compute new link lengths (They should be constant if solved correctly)
    for k=1:length(linkind) %loop over linkind
        del=xnode(cmat(linkind(k),1),i+1)-xnode(cmat(linkind(k),2),i+1);
        linksnew(k)=sqrt(real(del)^2+imag(del)^2); %compute length between nodes
    end
    linksnewtot=sum(linksnew);%Sums new link lengths (this should equal sum(links) during whole simulation)
    
    percdone=(i/tnum)*100.0; % Percentage of simulation calculated
    str=sprintf('Progress=%.2f %%: i=%i :Total Link Lengths=%.2f,',percdone,i,linksnewtot); %Make string with progress info
    disp(str)
end


%Section8: Animate Linkage
figsize=5.0; %relative figure size in inches 
xmin=min(min(real(xnode)));
xmax=max(max(real(xnode)));
ymin=min(min(imag(xnode)));
ymax=max(max(imag(xnode)));
hfig=figure(1); %initiate figure
figscale=figsize/max([xmax-xmin ymax-ymin]);
set(hfig,'units','inches','position',[2 2 (xmax-xmin)*figscale (ymax-ymin)*figscale]);
hax=axis(); %initiate axis
hold on
xlim([xmin xmax]) %Scale xaxis to fit max/min x of simulation
ylim([ymin ymax]) %Scale yaxis to fit max/min y of simulation
for k=1:tnum %loop over all timesteps
    hax = gca; %initialize current axis
    hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
    for i=1:length(cmat) %loop over links

            
            p1in=cmat(i,1); %link i node 1
            p2in=cmat(i,2); %link i node 2
            scatter(real([xnode(p1in,1:k),xnode(p2in,1:k)]),imag([xnode(p1in,1:k),xnode(p2in,1:k)]),'.k'); %Plot position of each node as black dot 
            
            %If groundlink, plot dashed line, else make solid
            if sum(i==gndlink)==1
                plot(real([xnode(p1in,k),xnode(p2in,k)]),imag([xnode(p1in,k),xnode(p2in,k)]),'--') %Plot dashed line between
            else
                plot(real([xnode(p1in,k),xnode(p2in,k)]),imag([xnode(p1in,k),xnode(p2in,k)])) %Plot line between
            end
            midpos=xnode(p1in,k)+0.5*(xnode(p2in,k)-xnode(p1in,k));
            text(real(midpos),imag(midpos),[ '$\raisebox{.5pt}{\bf{\raisebox{-.9pt} {' num2str(i) '}}}$'],'Interpreter','latex');

    end
    
    %Puts node labels at nodes
    for j=1:length(inodes)
        text(real(xnode(j,k)),imag(xnode(j,k)),[ '$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {' num2str(j) '}}}$'], 'Interpreter', 'latex')
    end
    
    hxlabel=xlabel('X position');
    hylabel=ylabel('Y Position');
    htitle=title('Mechanism Animation');
    formatplot(hfig,hax,hxlabel,hylabel,htitle)
    pause(0.002) %how long to pause plot in seconds
    if k~=tnum %If it reaches the end, it wont clear the plot
        cla %otherwise clear axis to plot next time step
    end
end


%Section9: Compute Velocity and Acceleration of nodes

%Velocity Calculation
velocity=zeros(length(inodes),tnum-1);
for i=1:tnum-1
    % complete this part yourself
    for j=1:length(inodes)
        if i==1
             velocity(j,i)=(xnode(j,i+1)-xnode(j,i))/dt;
        end
        %Use symmetric difference differentiation (less error)
        if i~=1
             velocity(j,i)=(xnode(j,i+1)-xnode(j,i-1))/(2*dt);
        end
       
    end
end


%Acceleration code goes here
acceleration=zeros(length(inodes),tnum-2);
for i=1:tnum-2
    %First order differentiation
    if i==1
        acceleration(:,i)=(velocity(:,i+1)-velocity(:,i))/dt;
    end
    %Use symmetric difference differentiation (less error)
    if i~=1
       acceleration(:,i)=(velocity(:,i+1)-velocity(:,i-1))/(2*dt);
    end
end

%Additional smoothening of the position/velocity data may need to be
%implemented if there are "spikes" in the acceleration. This has to do with
%the stability of the numerical procedure used?

 crankangles=zeros(1,tnum);
 crank_node=cmat(crlink,:);
 crankangles(1)=(180/pi*angle(inodes(crank_node(2))-inodes(crank_node(1))));
 if(crankangles(1)<0)
     crankangles(1)=crankangles(1)+360;
 end
 
 
 for j=2:(tnum)
     crankangles(j)=crankangles(j-1)+crangvel*dt;
      if(crankangles(j)<0)
        crankangles(j)=crankangles(j)+360;
      elseif(crankangles(j)>360)
        crankangles(j)=crankangles(j)-360;
      end
 end

times=zeros(1,tnum);

for j=2:tnum
    times(j)=times(j-1)+dt;
end


%position plot
pos_norm=zeros(1,tnum);
pos_norm(1:tnum)=real(xnode(1,:));




hfig=figure(2); %initiate figure
figscale=figsize/360;

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(pos_norm,imag(xnode(1,:)));

hxlabel=xlabel('Position X (in))');
hylabel=ylabel('Position Y (in)');
htitle=title('Position Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Position');


%position plot
pos_norm=zeros(1,tnum);
pos_norm(1:tnum)=real(xnode(1,:));




hfig=figure(11); %initiate figure
figscale=figsize/360;

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(crankangles,real(xnode(1,:)));

hxlabel=xlabel('Crank Angle(deg)');
hylabel=ylabel('Position X (in)');
htitle=title('Position X Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Position');

%position plot
pos_norm=zeros(1,tnum);
pos_norm(1:tnum)=real(xnode(1,:));




hfig=figure(10); %initiate figure
figscale=figsize/360;

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(crankangles,imag(xnode(1,:)));

hxlabel=xlabel('Crank Angle (deg))');
hylabel=ylabel('Position Y (in)');
htitle=title('Position Y Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Position');

%velocity y plot
velocity_norm=zeros(1,tnum);
velocity_norm(1:tnum-1)=imag(velocity(1,:));
velocity_norm(tnum)=velocity_norm(tnum-1);




hfig=figure(3); %initiate figure
figscale=figsize/360;

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(crankangles,velocity_norm);

hxlabel=xlabel('Crank Angle(deg)');
hylabel=ylabel('Velocity Y(in/s)');
htitle=title('Velocity Y Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Velocity');


%velocity y plot
velocity_norm=zeros(1,tnum);
velocity_norm(1:tnum-1)=real(velocity(1,:))*1/0.0254*5;
velocity_norm(tnum)=velocity_norm(tnum-1);




hfig=figure(4); %initiate figure
figscale=figsize/360;

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(crankangles,velocity_norm);

hxlabel=xlabel('Crank Angle(deg)');
hylabel=ylabel('Velocity X(fpm)');
htitle=title('Velocity X Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Velocity');



%acceleration plot
accel_norm=zeros(1,tnum);
accel_norm(1:tnum-2)=real(acceleration(1,:));
accel_norm(tnum-1)=accel_norm(tnum-2);
accel_norm(tnum)=accel_norm(tnum-1);



%acceleration x plot


hfig=figure(5); %initiate figure
figscale=20/360;
hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(crankangles,accel_norm);

hxlabel=xlabel('Crank Angles (deg)');
hylabel=ylabel('Acceleration X(in/s^2)');
htitle=title('Acceleration X Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Acceleration');





%acceleration y plot

accel_norm=zeros(1,tnum);
accel_norm(1:tnum-2)=imag(acceleration(1,:));
accel_norm(tnum-1)=accel_norm(tnum-2);
accel_norm(tnum)=accel_norm(tnum-1);

hfig=figure(6); %initiate figure
figscale=20/360;
hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(crankangles,accel_norm);

hxlabel=xlabel('Crank Angles(deg)');
hylabel=ylabel('Acceleration Y(in/s^2)');
htitle=title('Acceleration Y Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Acceleration');




%velocity plot
velocity_norm=zeros(1,tnum);
velocity_norm(1:tnum-1)=norm(velocity(1,:));
velocity_norm(tnum)=velocity_norm(tnum-1);



% velocity norm plot
hfig=figure(7); %initiate figure
figscale=figsize/360;

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(times,velocity_norm);

% hxlabel=xlabel('Time(s)');
% hylabel=ylabel('Velocity(m/s)');
% htitle=title('Velocity Plot');
% formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Velocity Norm');

% velocity(Vy vs Vx) plot
% velocity_norm=zeros(1,tnum);
% velocity_norm(1:tnum-1)=real(velocity(1,:));
% velocity_norm(tnum)=velocity_norm(tnum-1);
% 



hfig=figure(8); %initiate figure
figscale=figsize/360;

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(real(velocity(1,:)),imag(velocity(1,:)));

hxlabel=xlabel('Velocity(x)');
hylabel=ylabel('Velocity(y)');
htitle=title('Velocity Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Velocity');


% power plot

mass=0.7;
scale=0.0254;

power=zeros(1,tnum);
for(i=1:tnum-2)
    power(i)=abs(mass*real(velocity(1,i))*real(acceleration(1,i))*scale)+abs(mass*imag(velocity(1,i))*imag(acceleration(1,i))*scale)+abs(mass*9.81*imag(velocity(1,i)));
    power(i)=abs(power(i)*scale);
end
power(tnum-1)=mass*norm(velocity(1,tnum-1))*norm(acceleration(1,tnum-2))*scale^2;
power(tnum)=mass*norm(velocity(1,tnum-1))*norm(acceleration(1,tnum-2))*scale^2;

hfig=figure(9); %initiate figure
figscale=figsize/360;

hax=axis(); %initiate axis
hold on



hax = gca; %initialize current axis
hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
plot(crankangles,power);

hxlabel=xlabel('Angle');
hylabel=ylabel('Power');
htitle=title('Power Plot');
formatplot(hfig,hax,hxlabel,hylabel,htitle);
legend('Power');