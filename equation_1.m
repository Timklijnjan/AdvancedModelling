function [alphaddotcoeff_1,thetaddotcoeff_final_1,epsddotcoeff_1,Q_1,equation_1_final]=equation_1(alpha__hat,theta,epsilon,theta__dot,epsilon__dot,v,g)
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
w=v/r__a; %check whether it's true or not

I_w=1;
 x__dot= -v*sin(theta)*theta__dot;
% xddot= -v*(cos(theta)*(theta__dot)^2+sin(theta)*thetaddot);
 y__dot= v*cos(theta)*theta__dot;
 %yddot=  v*(-sin(theta)*(theta__dot)^2+cos(theta)*thetaddot);
%alphaddotcoeff_1 contains as variables only s__1,s__2
F_alphaddotcoeff_1=@(s__1,s__2) (r__a^2 - 2*r__a*s__2 + s__1^2 + s__2^2)*M__a + 2*I__yy;
alphaddotcoeff_1=F_alphaddotcoeff_1(s__1,s__2);

%thetaddotcoeff_1 contains as variables only alpha_hat, epsilon, s__1,s__2
F_thetaddotcoeff_1=@(alpha__hat,epsilon,s__1,s__2) -r__a*(-s__2 + r__a)*sin(epsilon)*M__a*cos(alpha__hat) + M__a*r__a*s__1*sin(epsilon)*sin(alpha__hat) + M__a*(r__a^2 - 2*r__a*s__2 + s__1^2 + s__2^2)*sin(epsilon) + 2*I__yz;
thetaddotcoeff_1=F_thetaddotcoeff_1(alpha__hat,epsilon,s__1,s__2);

%epsddotcoeff is constant
epsddotcoeff_1=2*I__xy;

%xddotcoeff_1 contains as variables only alpha_hat, epsilon, s__1,s__2
F_xddotcoeff_1=@(alpha__hat,theta,epsilon,s__1,s__2) (-s__1*sin(theta)*sin(epsilon) + (-s__2 + r__a)*cos(theta))*M__a*cos(alpha__hat) + M__a*(-(-s__2 + r__a)*sin(theta)*sin(epsilon) - s__1*cos(theta))*sin(alpha__hat);
xddotcoeff_1=F_xddotcoeff_1(alpha__hat,theta,epsilon,s__1,s__2);

%thetaddotextra_x_1 contains as variables only theta and xddotcoeff_1, so
%also variables for xddotcoeff_1 have to be put
F_thetaddotextra_x_1=@(alpha__hat,theta,epsilon,s__1,s__2) -v*sin(theta)*F_xddotcoeff_1(alpha__hat,theta,epsilon,s__1,s__2);
thetaddotextra_x_1=F_thetaddotextra_x_1(alpha__hat,theta,epsilon,s__1,s__2);

%yddotcoeff_1 contains as variables only alpha_hat,theta, epsilon, s__1,s__2
F_yddotcoeff_1=@(alpha_hat,theta, epsilon, s__1,s__2) (s__1*cos(theta)*sin(epsilon) + (-s__2 + r__a)*sin(theta))*M__a*cos(alpha__hat) + M__a*((-s__2 + r__a)*cos(theta)*sin(epsilon) - s__1*sin(theta))*sin(alpha__hat);
yddotcoeff_1=F_yddotcoeff_1(alpha__hat,theta, epsilon, s__1,s__2);

%thetaddotextra_x_1 contains as variables only theta and yddotcoeff_1, so
%also variables for yddotcoeff_1 have to be put
F_thetaddotextra_y_1=@(alpha__hat,theta, epsilon, s__1,s__2) v*cos(theta)*F_yddotcoeff_1(alpha__hat,theta, epsilon, s__1,s__2);
thetaddotextra_y_1=F_thetaddotextra_y_1(alpha__hat,theta, epsilon, s__1,s__2);

%thetaddotcoeff_final_1 contains as variables only thetaddotcoeff_1,thetaddotextra_x_1,thetaddotextra_y_1;
F_thetaddotcoeff_final_1=@(alpha__hat,theta, epsilon, s__1,s__2) F_thetaddotcoeff_1(alpha__hat,epsilon,s__1,s__2) +F_thetaddotextra_x_1(alpha__hat,theta,epsilon,s__1,s__2) +F_thetaddotextra_y_1(alpha__hat,theta, epsilon, s__1,s__2);
thetaddotcoeff_final_1=F_thetaddotcoeff_final_1(alpha__hat,theta, epsilon, s__1,s__2);% still to do+s__1*y__ddot*cos(theta)*sin(epsilon);

%equation_1 contains as variables alpha__hat,theta, epsilon, s__1,s__2,theta__dot,epsilon__dot
F_equation_1=@(alpha__hat,theta, epsilon, s__1,s__2,theta__dot,epsilon__dot) -2*(theta__dot*(-s__2 + r__a)*cos(epsilon) + epsilon__dot*s__1)*(theta__dot*s__1*cos(epsilon) - epsilon__dot*(-s__2 + r__a))*M__a*cos(alpha__hat)^2 + (-(theta__dot*(r__a + s__1 - s__2)*cos(epsilon) - epsilon__dot*(r__a - s__1 - s__2))*(theta__dot*(r__a - s__1 - s__2)*cos(epsilon) + epsilon__dot*(r__a + s__1 - s__2))*sin(alpha__hat) + cos(epsilon)^2*r__a*s__1*theta__dot^2 + (-2*r__a^2*epsilon__dot*theta__dot + 2*r__a*s__2*epsilon__dot*theta__dot + g*s__1)*cos(epsilon) - r__a*s__1*(epsilon__dot^2 + theta__dot^2))*M__a*cos(alpha__hat) + M__a*(theta__dot^2*(-s__2 + r__a)*r__a*cos(epsilon)^2 + ((2*s__1*epsilon__dot*theta__dot + g)*r__a - g*s__2)*cos(epsilon) - r__a*(epsilon__dot^2 + theta__dot^2)*(-s__2 + r__a))*sin(alpha__hat) + theta__dot^2*M__a*s__1*(-s__2 + r__a)*cos(epsilon)^2 + 2*M__a*cos(epsilon)*s__1^2*epsilon__dot*theta__dot + (-r__a*s__1*epsilon__dot^2 + s__1*s__2*epsilon__dot^2)*M__a;
equation_1=F_equation_1(alpha__hat,theta, epsilon, s__1,s__2,theta__dot,epsilon__dot);

%equationextra_1 contains as variables theta, xddotcoeff,yddotcoeff
F_equationextra_1=@(alpha__hat,theta,epsilon,s__1,s__2)-v*cos(theta)*(theta__dot)^2*F_xddotcoeff_1(alpha__hat,theta,epsilon,s__1,s__2) -v*sin(theta)*(theta__dot)^2*F_yddotcoeff_1(alpha__hat,theta,epsilon,s__1,s__2);
equationextra_1=F_equationextra_1(alpha__hat,theta,epsilon,s__1,s__2);

F_equation_1_final=@(alpha__hat,theta, epsilon, s__1,s__2,theta__dot,epsilon__dot) F_equation_1(alpha__hat,theta, epsilon, s__1,s__2,theta__dot,epsilon__dot)+F_equationextra_1(alpha__hat,theta,epsilon,s__1,s__2);
equation_1_final=F_equation_1_final(alpha__hat,theta, epsilon, s__1,s__2,theta__dot,epsilon__dot);

%Q_1 contains as variables only theta,epsilon,epsilon__dot;
F_Q_1=@(theta,epsilon,epsilon__dot) I_w*w*(cos(theta)*cos(epsilon)*theta - sin(theta)*sin(epsilon)*epsilon__dot);
%Q_1=F_Q_1(theta,epsilon,epsilon__dot);
Q_1=0;
end