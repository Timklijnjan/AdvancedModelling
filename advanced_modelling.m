%-----------------------------------------------------------------------------------------------------
% PROBLEM FORMULATION
% The equation i is in the form:
% alphaddot*(alphaddotcoeff_i)+thetaddot*(thetaddotcoeff_final_i)+epsddot*(epsddotcoeff_i)=Q_i-equation_i_final
% The system is then:
% A*uddot=f 
% uddot is 3x1 vector [alphaddot;thetaddot; epsilonddot]
% A is 3x3 matrix, containing udot, u and constants
% f is 3x1 rhs, containing udot and u 
%-----------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------
% SOLVING STRATEGY
% system will be solved using Euler-forward, with timestep dt and inital
% conditions: alpha__hat(0), theta(0), epsilon(0), alpha__hatdot(0),
% theta__dot(0), epsilon__dot(0); i.e. u(0) and udot(0)
% This approach leads to solving the system time-step-wise, i.e.
% uddot(0)=A(0)\f(0);
% udot(1)=udot(0)+dt*uddot(0);
% u(1)= u(0)+ dt*udot(0);
% uddot(1)= A(1)\f(1); and so on
%-----------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------
% PARAMETERS
% I_w moment of inertia about the axis going normal through the center of the wheel
%  w angular velocity of the wheel
%  r__a radius of the wheel
%  s__1 inital position of the centre of mass on x-direction
%  s__2 inital position of the centre of mass on y-direction
%  M__a mass of the bike
%  I__xx,I__xy,I__xz,I__yy, I__yz, I__zz inertial momentum 
%  v velocity of the bike
%  x__dot=v*cos(theta), y__dot=v*sin(theta)
%-----------------------------------------------------------------------------------------------------
%PARAMETERS
r__a=0.32; %Rtilde, check whether is the radius of wheel or something else
s__1=0;
s__2=1;
M__a=81; %m1+m2, check whether its m1 / m2/ m1+m2
I__xx=60.57;
I__xy=0;
I__xz=0;
I__yy=60.62;
I__yz=0;
I__zz=0.15;
I_w=0.09;
v=1;
w=v/r__a; %check whether it's true or not
g=9.81;

%----------------------------------------------------------------------------------------------------------------------------%
%RESOLUTION
n=1000;
dt=0.01;
x=zeros(1,n+1);
y=zeros(1,n+1);
udot=zeros(3,n+1);
uddot=zeros(3,n+1);
u=zeros(3,n+1);
u(:,1)=[0;0;0]; %initial condition (alpha(0),theta(0),epsilon(0))
udot(:,1)=[0;pi/5;0]; %initial condition (alphadot(0),thetadot(0),epsilondot(0))
for i=1:n
    [alphaddotcoeff_1,thetaddotcoeff_final_1,epsddotcoeff_1,Q_1,equation_1_final]=equation_1(u(1,i),u(2,i),u(3,i),udot(2,i),udot(3,i),v,g);
    [alphaddotcoeff_2,thetaddotcoeff_final_2,epsddotcoeff_2,Q_2,equation_2_final]=equation_2(u(1,i),u(2,i),u(3,i),udot(1,i),udot(2,i),udot(3,i),v,g);
    [alphaddotcoeff_3,thetaddotcoeff_final_3,epsddotcoeff_3,Q_3,equation_3_final]=equation_3(u(1,i),u(2,i),u(3,i),udot(1,i),udot(2,i),udot(3,i),v,g);
    A=[alphaddotcoeff_1 thetaddotcoeff_final_1 epsddotcoeff_1; alphaddotcoeff_2 thetaddotcoeff_final_2 epsddotcoeff_2; alphaddotcoeff_3 thetaddotcoeff_final_3 epsddotcoeff_3];
    b=[Q_1-equation_1_final;Q_2-equation_2_final;Q_3-equation_3_final];

uddot(:,i)=A\b;
udot(:,i+1)=udot(:,i)+dt*uddot(:,i);
u(:,i+1)=u(:,i)+dt*udot(:,i+1);
x(1,i+1)=x(1,i)+dt*v*cos(u(2,i));
y(1,i+1)=y(1,i)+dt*v*sin(u(2,i));
end


%------------------------------------------------------------------------------------------------------------------------------%
%PLOTS
figure(1);

subplot(2,2,1);
plot(0:dt:n*dt,rad2deg(u(1,:)),'o-')
% Ax = gca;
% Ax.TickLabelInterpreter = 'latex';
% yt = Ax.YTick;
% ytl = {'$-\pi/2$', '$-\pi/4$', '0', '$\pi/4$', '$\pi/2$'};
% ytv = linspace(min(yt), max(yt), numel(ytl));
% set(Ax, 'YTick',ytv, 'YTickLabel',ytl)
xlabel('$t$', 'Interpreter','latex');
ylabel('$\alpha$', 'Interpreter','latex');

subplot(2,2,2);
plot(0:dt:n*dt,rad2deg(u(2,:)),'o-')
% Ax = gca;
% Ax.TickLabelInterpreter = 'latex';
% yt = Ax.YTick;
% ytl = {'$-\pi/2$', '$-\pi/4$', '0', '$\pi/4$', '$\pi/2$'};
% ytv = linspace(min(yt), max(yt), numel(ytl));
% set(Ax, 'YTick',ytv, 'YTickLabel',ytl)
xlabel('$t$', 'Interpreter','latex');
ylabel('$\theta$', 'Interpreter','latex');

subplot(2,2,3);
plot(0:dt:n*dt,rad2deg(u(3,:)),'o-')
% Ax = gca;
% Ax.TickLabelInterpreter = 'latex';
% yt = Ax.YTick;
% ytl = {'$-\pi/2$', '$-\pi/4$', '0', '$\pi/4$', '$\pi/2$'};
% ytv = linspace(min(yt), max(yt), numel(ytl));
% set(Ax, 'YTick',ytv, 'YTickLabel',ytl)
xlabel('$t$', 'Interpreter','latex');
ylabel('$\epsilon$', 'Interpreter','latex');

figure(2);
xlocs ='bike';
Unicyclemoviemaker(x,y,u(2,:),u(3,:),u(1,:),r__a,s__2,xlocs)
