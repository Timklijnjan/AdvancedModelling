%-----------------------------------------------------------------------------------------------------
% PROBLEM FORMULATION
% The equation i is in the form:
% alphaddot*(alphaddotcoeff_i)+thetaddot*(thetaddotcoeff_final_i)+epsddot*(epsddotcoeff_i)=Q_i-equation_i_final
% The system is then:
% A*uddot=f 
% uddot is 3x1 vector [alphaddot;thetaddot; epsilonddot]
% A is 3x3 matrix, containing udot and u
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
%  x=v*cos(theta), y=v*sin(theta)
%-----------------------------------------------------------------------------------------------------

x__dot=-v*sin(theta)*theta__dot;
xddot=-v*(cos(theta)*(theta__dot)^2+sin(theta)*thetaddot);
y__dot=v*cos(theta)*theta__dot;
yddot=v*(-sin(theta)*(theta__dot)^2+cos(theta)*thetaddot);

%equation with derivatives w.r.t alphaddot --> 1
alphaddotcoeff_1=(r__a^2 - 2*r__a*s__2 + s__1^2 + s__2^2)*M__a + 2*I__yy;
thetaddotcoeff_1=-r__a*(-s__2 + r__a)*sin(epsilon)*M__a*cos(alpha__hat) + M__a*r__a*s__1*sin(epsilon)*sin(alpha__hat) + M__a*(r__a^2 - 2*r__a*s__2 + s__1^2 + s__2^2)*sin(epsilon) + 2*I__yz;
epsddotcoeff_1=2*I__xy;
xddotcoeff_1=(-s__1*sin(theta)*sin(epsilon) + (-s__2 + r__a)*cos(theta))*M__a*cos(alpha__hat) + M__a*(-(-s__2 + r__a)*sin(theta)*sin(epsilon) - s__1*cos(theta))*sin(alpha__hat);
thetaddotextra_x_1=-v*sin(theta)*xddotcoeff_1;
yddotcoeff_1=(s__1*cos(theta)*sin(epsilon) + (-s__2 + r__a)*sin(theta))*M__a*cos(alpha__hat) + M__a*((-s__2 + r__a)*cos(theta)*sin(epsilon) - s__1*sin(theta))*sin(alpha__hat);
thetaddotextra_y_1=v*cos(theta)*yddotcoeff_1;
thetaddotcoeff_final_1=thetaddotcoeff_1+thetaddotextra_x_1+thetaddotextra_y_1;
equation_1=-2*(theta__dot*(-s__2 + r__a)*cos(epsilon) + epsilon__dot*s__1)*(theta__dot*s__1*cos(epsilon) - epsilon__dot*(-s__2 + r__a))*M__a*cos(alpha__hat)^2 + (-(theta__dot*(r__a + s__1 - s__2)*cos(epsilon) - epsilon__dot*(r__a - s__1 - s__2))*(theta__dot*(r__a - s__1 - s__2)*cos(epsilon) + epsilon__dot*(r__a + s__1 - s__2))*sin(alpha__hat) + cos(epsilon)^2*r__a*s__1*theta__dot^2 + (-2*r__a^2*epsilon__dot*theta__dot + 2*r__a*s__2*epsilon__dot*theta__dot + g*s__1)*cos(epsilon) + (s__1*y__ddot*cos(theta) - s__1*x__ddot *sin(theta))*sin(epsilon) + x__ddot*(-s__2 + r__a)*cos(theta) + y__ddot*(-s__2 + r__a)*sin(theta) - r__a*s__1*(epsilon__dot^2 + theta__dot^2))*M__a*cos(alpha__hat) + M__a*(theta__dot^2*(-s__2 + r__a)*r__a*cos(epsilon)^2 + ((2*s__1*epsilon__dot*theta__dot + g)*r__a - g*s__2)*cos(epsilon) + (y__ddot*(-s__2 + r__a)*cos(theta) - x__ddot*(-s__2 + r__a)*sin(theta))*sin(epsilon) - s__1*x__ddot *cos(theta) - s__1*y__ddot*sin(theta) - r__a*(epsilon__dot^2 + theta__dot^2)*(-s__2 + r__a))*sin(alpha__hat) + theta__dot^2*M__a*s__1*(-s__2 + r__a)*cos(epsilon)^2 + 2*M__a*cos(epsilon)*s__1^2*epsilon__dot*theta__dot + (-r__a*s__1*epsilon__dot^2 + s__1*s__2*epsilon__dot^2)*M__a;
equationextra_1=-v*cos(theta)*(theta__dot)^2*xddotcoeff_1 -v*sin(theta)*(theta_dot)^2*yddotcoeff_1;
equation_1_final=equation_1+equationextra_1;
Q_1=I_w*w(cos(theta)*cos(epsilon)*theta - sin(theta)*sin(epsilon)*epsilon__dot);

%equation with derivatives w.r.t thetaddot --> 2
alphaddotcoeff_2=-(-s__2 + r__a)*r__a*sin(epsilon)*M__a*cos(alpha__hat) - M__a*(-r__a*sin(alpha__hat)*s__1 - r__a^2 + 2*s__2*r__a - s__1^2 - s__2^2)*sin(epsilon) + 2*I__yz;
thetaddotcoeff_2=2*(-1/2*r__a^2 + s__2*r__a + 1/2*s__1^2 - 1/2*s__2^2)*cos(epsilon)^2*M__a*cos(alpha__hat)^2 + (((2*r__a*s__1 - 2*s__1*s__2)*sin(alpha__hat) - 2*r__a*(-r__a + s__2))*cos(epsilon)^2 + 2*r__a*(-r__a + s__2))*M__a*cos(alpha__hat) - 2*(r__a*sin(alpha__hat)*s__1 + r__a^2/2 + s__1^2/2)*M__a*cos(epsilon)^2 + 2*M__a*r__a*s__1*sin(alpha__hat) + 2*(r__a^2 - s__2*r__a + 1/2*s__1^2 + 1/2*s__2^2)*M__a + 2*I__zz;
epsddotcoeff_2=2*(-r__a*s__1 + s__1*s__2)*cos(epsilon)*M__a*cos(alpha__hat)^2 + ((-r__a^2 + 2*r__a*s__2 + s__1^2 - s__2^2)*sin(alpha__hat) + r__a*s__1)*cos(epsilon)*M__a*cos(alpha__hat) + 2*M__a*(r__a*sin(alpha__hat)/2 + s__1/2)*(-s__2 + r__a)*cos(epsilon) + 2*I__xz;
xddotcoeff_2=((-s__2 + r__a)*cos(theta)*sin(epsilon) - s__1*sin(theta))*M__a*cos(alpha__hat) - M__a*(s__1*cos(theta)*sin(alpha__hat) + cos(theta)*r__a)*sin(epsilon) - M__a*(-s__2 + r__a)*sin(theta)*sin(alpha__hat);
yddotcoeff_2=((-s__2 + r__a)*sin(theta)*sin(epsilon) + s__1*cos(theta))*M__a*cos(alpha__hat) - M__a*(s__1*sin(theta)*sin(alpha__hat) + sin(theta)*r__a)*sin(epsilon) + M__a*(-s__2 + r__a)*cos(theta)*sin(alpha__hat);
thetaddotextra_x_2=-v*sin(theta)*xddotcoeff_2;
thetaddotextra_y_2=v*cos(theta)*yddotcoeff_2;
thetaddotcoeff_final_2=thetaddotcoeff_2+thetaddotextra_x_2+thetaddotextra_y_2;
equation_2=2*((2*r__a*s__1*alpha__hatdot*theta__dot - 2*s__1*s__2*alpha__hatdot*theta__dot)*cos(epsilon)^2 + (theta__dot*epsilon__dot*(r__a + s__1 - s__2)*(r__a - s__1 - s__2)*sin(epsilon) - epsilon__dot*alpha__hatdot*r__a^2 + 2*epsilon__dot*alpha__hatdot*r__a*s__2 + epsilon__dot*alpha__hatdot*s__1^2 - epsilon__dot*alpha__hatdot*s__2^2)*cos(epsilon) + epsilon__dot^2*sin(epsilon)*s__1*(-s__2 + r__a))*M__a*cos(alpha__hat)^2 + (((2*r__a^2*alpha__hatdot*theta__dot - 4*r__a*s__2*alpha__hatdot*theta__dot - 2*s__1^2*alpha__hatdot*theta__dot + 2*s__2^2*alpha__hatdot*theta__dot)*sin(alpha__hat) - 2*r__a*s__1*alpha__hatdot*theta__dot)*cos(epsilon)^2 + (-4*theta__dot*epsilon__dot*(-s__2 + r__a)*(s__1*sin(alpha__hat) + r__a)*sin(epsilon) + (4*r__a*s__1*alpha__hatdot*epsilon__dot - 4*s__1*s__2*alpha__hatdot*epsilon__dot)*sin(alpha__hat))*cos(epsilon) + (epsilon__dot^2*(r__a - s__1 - s__2)*(r__a + s__1 - s__2)*sin(alpha__hat) + `x__ddot `*(-s__2 + r__a)*cos(theta) + `y__ddot `*(-s__2 + r__a)*sin(theta) - r__a*(-alpha__hatdot^2 + epsilon__dot^2)*s__1)*sin(epsilon) + s__1*`y__ddot `*cos(theta) - s__1*`x__ddot `*sin(theta) + 2*r__a*s__1*alpha__hatdot*theta__dot)*M__a*cos(alpha__hat) - 2*(r__a*(r__a*alpha__hatdot*theta__dot - s__2*alpha__hatdot*theta__dot)*sin(alpha__hat) + r__a*s__1*alpha__hatdot*theta__dot - s__1*alpha__hatdot*theta__dot*s__2)*M__a*cos(epsilon)^2 + 2*M__a*(theta__dot*epsilon__dot*(2*r__a*sin(alpha__hat)*s__1 + r__a^2 + s__1^2)*sin(epsilon) + (r__a*alpha__hatdot*epsilon__dot - s__2*alpha__hatdot*epsilon__dot)*(-s__2 + r__a))*cos(epsilon) - M__a*((s__1*`x__ddot `*cos(theta) + s__1*`y__ddot `*sin(theta) - ((alpha__hatdot^2 - epsilon__dot^2)*r__a + s__2*(-alpha__hatdot^2 + epsilon__dot^2))*r__a)*sin(alpha__hat) + `x__ddot `*cos(theta)*r__a + `y__ddot `*sin(theta)*r__a + epsilon__dot^2*s__1*r__a - s__1*s__2*epsilon__dot^2)*sin(epsilon) - M__a*(-`y__ddot `*(-s__2 + r__a)*cos(theta) + `x__ddot `*(-s__2 + r__a)*sin(theta) - 2*r__a*(r__a*alpha__hatdot*theta__dot - s__2*alpha__hatdot*theta__dot))*sin(alpha__hat);
equationextra_2=-v*cos(theta)*(theta__dot)^2*xddotcoeff_2 -v*sin(theta)*(theta_dot)^2*yddotcoeff_2;
equation_2_final=equation_2+equationextra_2;
Q_2=I_w*w(sin(theta)*cos(epsilon)*theta + cos(theta)*sin(epsilon)*epsilon__dot);

%equation with derivatives w.r.t epsilonddott --> 3
alphaddotcoeff_3=2*I__xy;
thetaddotcoeff_3=-(2*r__a*s__1 - 2*s__1*s__2)*cos(epsilon)*M__a*cos(alpha__hat)^2 + ((-r__a^2 + 2*r__a*s__2 + s__1^2 - s__2^2)*sin(alpha__hat) + r__a*s__1)*cos(epsilon)*M__a*cos(alpha__hat) - M__a*((-r__a^2 + r__a*s__2)*sin(alpha__hat) - r__a*s__1 + s__1*s__2)*cos(epsilon) + 2*I__xz;
epsddotcoeff_3=-(-r__a^2 + 2*r__a*s__2 + s__1^2 - s__2^2)*M__a*cos(alpha__hat)^2 + ((-2*r__a*s__1 + 2*s__1*s__2)*sin(alpha__hat) + 2*r__a*(s__2 - r__a))*M__a*cos(alpha__hat) + 2*r__a*s__1*M__a*sin(alpha__hat) + (r__a^2 + s__1^2)*M__a + 2*I__xx;
xddotcoeff_3=(-s__2 + r__a)*sin(theta)*cos(epsilon)*M__a*cos(alpha__hat) - M__a*(sin(theta)*s__1*sin(alpha__hat) + sin(theta)*r__a)*cos(epsilon);
yddotcoeff_3=-(-s__2 + r__a)*cos(theta)*cos(epsilon)*M__a*cos(alpha__hat) - M__a*(-cos(theta)*s__1*sin(alpha__hat) - cos(theta)*r__a)*cos(epsilon);
thetaddotextra_x_3=-v*sin(theta)*xddotcoeff_3;
thetaddotextra_y_3=v*cos(theta)*yddotcoeff_3;
thetaddotcoeff_final_3=thetaddotcoeff_3+thetaddotextra_x_3+thetaddotextra_y_3;
equation_3=-((theta__dot^2*(r__a + s__1 - s__2)*(r__a - s__1 - s__2)*sin(epsilon) + 2*r__a^2*alpha__hatdot*theta__dot - 4*r__a*alpha__hatdot*theta__dot*s__2 - 2*alpha__hatdot*theta__dot*s__1^2 + 2*alpha__hatdot*theta__dot*s__2^2)*cos(epsilon) + 4*r__a*s__1*alpha__hatdot*epsilon__dot - 4*s__1*alpha__hatdot*epsilon__dot*s__2)*M__a*cos(alpha__hat)^2 + (((2*theta__dot^2*s__1*(-s__2 + r__a)*sin(epsilon) + 4*r__a*s__1*alpha__hatdot*theta__dot - 4*s__1*alpha__hatdot*theta__dot*s__2)*sin(alpha__hat) + 2*theta__dot^2*r__a*(-s__2 + r__a)*sin(epsilon) - `y__ddot `*(-s__2 + r__a)*cos(theta) + `x__ddot `*(-s__2 + r__a)*sin(theta) + 2*r__a*(r__a*alpha__hatdot*theta__dot - s__2*alpha__hatdot*theta__dot))*cos(epsilon) + (-2*r__a^2*alpha__hatdot*epsilon__dot + 4*r__a*s__2*alpha__hatdot*epsilon__dot + 2*s__1^2*alpha__hatdot*epsilon__dot - 2*s__2^2*alpha__hatdot*epsilon__dot)*sin(alpha__hat) + g*(-s__2 + r__a)*sin(epsilon) + 2*r__a*s__1*alpha__hatdot*epsilon__dot)*M__a*cos(alpha__hat) - M__a*((2*theta__dot^2*sin(epsilon)*r__a*s__1 + 2*r__a*s__1*alpha__hatdot*theta__dot - s__1*`y__ddot `*cos(theta) + s__1*`x__ddot `*sin(theta))*sin(alpha__hat) + theta__dot^2*(r__a^2 + s__1^2)*sin(epsilon) + 2*alpha__hatdot*theta__dot*s__1^2 + sin(theta)*`x__ddot `*r__a - cos(theta)*`y__ddot `*r__a)*cos(epsilon) - (g*sin(epsilon)*s__1 - 2*r__a*(r__a*alpha__hatdot*epsilon__dot - s__2*alpha__hatdot*epsilon__dot))*M__a*sin(alpha__hat) - M__a*sin(epsilon)*g*r__a + (2*r__a*s__1*alpha__hatdot*epsilon__dot - 2*s__1*s__2*alpha__hatdot*epsilon__dot)*M__a;
equationextra_3=-v*cos(theta)*(theta__dot)^2*xddotcoeff_3 -v*sin(theta)*(theta_dot)^2*yddotcoeff_3;
equation_3_final=equation_3+equationextra_3;
Q_3= I_w*w*epsilon__dot*cos(epsilon);

%formulation of the system
A=[alphaddotcoeff_1 thetaddotcoeff_final_1 epsddotcoeff_1; alphaddotcoeff_2 thetaddotcoeff_final_2 epsddotcoeff_2; alphaddotcoeff_3 thetaddotcoeff_final_3 epsddotcoeff_3];
b=[Q_1-equationextra_1;Q_2-equationextra_2;Q_3-equationextra_3];
