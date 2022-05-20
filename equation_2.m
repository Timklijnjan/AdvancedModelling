function [alphaddotcoeff_2,thetaddotcoeff_final_2,epsddotcoeff_2,Q_2,equation_2_final]=equation_2(alpha__hat,theta,epsilon,alpha__hatdot,theta__dot,epsilon__dot,v,g)
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
%alphaddotcoeff_1 contains as variables only alpha__hat,epsilon,s__1,s__2
F_alphaddotcoeff_2=@(alpha__hat,epsilon,s__1,s__2) -(-s__2 + r__a)*r__a*sin(epsilon)*M__a*cos(alpha__hat) - M__a*(-r__a*sin(alpha__hat)*s__1 - r__a^2 + 2*s__2*r__a - s__1^2 - s__2^2)*sin(epsilon) + 2*I__yz;
alphaddotcoeff_2=F_alphaddotcoeff_2(alpha__hat,epsilon,s__1,s__2);

%thetaddotcoeff_2 contains as variables only alpha__hat,epsilon,s__1,s__2
F_thetaddotcoeff_2=@(alpha__hat,epsilon,s__1,s__2) 2*(-1/2*r__a^2 + s__2*r__a + 1/2*s__1^2 - 1/2*s__2^2)*cos(epsilon)^2*M__a*cos(alpha__hat)^2 + (((2*r__a*s__1 - 2*s__1*s__2)*sin(alpha__hat) - 2*r__a*(-r__a + s__2))*cos(epsilon)^2 + 2*r__a*(-r__a + s__2))*M__a*cos(alpha__hat) - 2*(r__a*sin(alpha__hat)*s__1 + r__a^2/2 + s__1^2/2)*M__a*cos(epsilon)^2 + 2*M__a*r__a*s__1*sin(alpha__hat) + 2*(r__a^2 - s__2*r__a + 1/2*s__1^2 + 1/2*s__2^2)*M__a + 2*I__zz;
thetaddotcoeff_2=F_thetaddotcoeff_2(alpha__hat,epsilon,s__1,s__2);

%epsddotcoeff_2 contains as variables only alpha__hat,epsilon,s__1,s__2
F_epsddotcoeff_2=@(alpha__hat,epsilon,s__1,s__2) 2*(-r__a*s__1 + s__1*s__2)*cos(epsilon)*M__a*cos(alpha__hat)^2 + ((-r__a^2 + 2*r__a*s__2 + s__1^2 - s__2^2)*sin(alpha__hat) + r__a*s__1)*cos(epsilon)*M__a*cos(alpha__hat) + 2*M__a*(r__a*sin(alpha__hat)/2 + s__1/2)*(-s__2 + r__a)*cos(epsilon) + 2*I__xz;
epsddotcoeff_2=F_epsddotcoeff_2(alpha__hat,epsilon,s__1,s__2);

%xddotcoeff_2 contains as variables only alpha__hat,,theta,epsilon,s__1,s__2
F_xddotcoeff_2=@(alpha__hat,theta,epsilon,s__1,s__2) ((-s__2 + r__a)*cos(theta)*sin(epsilon) - s__1*sin(theta))*M__a*cos(alpha__hat) - M__a*(s__1*cos(theta)*sin(alpha__hat) + cos(theta)*r__a)*sin(epsilon) - M__a*(-s__2 + r__a)*sin(theta)*sin(alpha__hat);
xddotcoeff_2=F_xddotcoeff_2(alpha__hat,theta,epsilon,s__1,s__2);

%yddotcoeff_2 contains as variables only alpha__hat,,theta,epsilon,s__1,s__2
F_yddotcoeff_2=@(alpha__hat,theta,epsilon,s__1,s__2) ((-s__2 + r__a)*sin(theta)*sin(epsilon) + s__1*cos(theta))*M__a*cos(alpha__hat) - M__a*(s__1*sin(theta)*sin(alpha__hat) + sin(theta)*r__a)*sin(epsilon) + M__a*(-s__2 + r__a)*cos(theta)*sin(alpha__hat);
yddotcoeff_2=F_yddotcoeff_2(alpha__hat,theta,epsilon,s__1,s__2);

%thetaddotcoeff_2 contains as variables only alpha__hat,,theta,epsilon,s__1,s__2
F_thetaddotextra_x_2=@(alpha__hat,theta,epsilon,s__1,s__2) -v*sin(theta)*F_xddotcoeff_2(alpha__hat,theta,epsilon,s__1,s__2);
thetaddotextra_x_2=F_thetaddotextra_x_2(alpha__hat,theta,epsilon,s__1,s__2);

F_thetaddotextra_y_2=@(alpha__hat,theta,epsilon,s__1,s__2) v*cos(theta)*F_yddotcoeff_2(alpha__hat,theta,epsilon,s__1,s__2);
thetaddotextra_y_2=F_thetaddotextra_y_2(alpha__hat,theta,epsilon,s__1,s__2);

F_thetaddotcoeff_final_2=@(alpha__hat,theta,epsilon,s__1,s__2) F_thetaddotcoeff_2(alpha__hat,epsilon,s__1,s__2)+F_thetaddotextra_x_2(alpha__hat,theta,epsilon,s__1,s__2)+F_thetaddotextra_y_2(alpha__hat,theta,epsilon,s__1,s__2);
thetaddotcoeff_final_2=F_thetaddotcoeff_final_2(alpha__hat,theta,epsilon,s__1,s__2);

F_equation_2=@(alpha__hat,theta,epsilon,alpha__hatdot,theta__dot,epsilon__dot,s__1,s__2) 2*((2*r__a*s__1*alpha__hatdot*theta__dot - 2*s__1*s__2*alpha__hatdot*theta__dot)*cos(epsilon)^2 + (theta__dot*epsilon__dot*(r__a + s__1 - s__2)*(r__a - s__1 - s__2)*sin(epsilon) - epsilon__dot*alpha__hatdot*r__a^2 + 2*epsilon__dot*alpha__hatdot*r__a*s__2 + epsilon__dot*alpha__hatdot*s__1^2 - epsilon__dot*alpha__hatdot*s__2^2)*cos(epsilon) + epsilon__dot^2*sin(epsilon)*s__1*(-s__2 + r__a))*M__a*cos(alpha__hat)^2 + (((2*r__a^2*alpha__hatdot*theta__dot - 4*r__a*s__2*alpha__hatdot*theta__dot - 2*s__1^2*alpha__hatdot*theta__dot + 2*s__2^2*alpha__hatdot*theta__dot)*sin(alpha__hat) - 2*r__a*s__1*alpha__hatdot*theta__dot)*cos(epsilon)^2 + (-4*theta__dot*epsilon__dot*(-s__2 + r__a)*(s__1*sin(alpha__hat) + r__a)*sin(epsilon) + (4*r__a*s__1*alpha__hatdot*epsilon__dot - 4*s__1*s__2*alpha__hatdot*epsilon__dot)*sin(alpha__hat))*cos(epsilon) + (epsilon__dot^2*(r__a - s__1 - s__2)*(r__a + s__1 - s__2)*sin(alpha__hat) +  - r__a*(-alpha__hatdot^2 + epsilon__dot^2)*s__1)*sin(epsilon) + 2*r__a*s__1*alpha__hatdot*theta__dot)*M__a*cos(alpha__hat) - 2*(r__a*(r__a*alpha__hatdot*theta__dot - s__2*alpha__hatdot*theta__dot)*sin(alpha__hat) + r__a*s__1*alpha__hatdot*theta__dot - s__1*alpha__hatdot*theta__dot*s__2)*M__a*cos(epsilon)^2 + 2*M__a*(theta__dot*epsilon__dot*(2*r__a*sin(alpha__hat)*s__1 + r__a^2 + s__1^2)*sin(epsilon) + (r__a*alpha__hatdot*epsilon__dot - s__2*alpha__hatdot*epsilon__dot)*(-s__2 + r__a))*cos(epsilon) - M__a*(( - ((alpha__hatdot^2 - epsilon__dot^2)*r__a + s__2*(-alpha__hatdot^2 + epsilon__dot^2))*r__a)*sin(alpha__hat) + epsilon__dot^2*s__1*r__a - s__1*s__2*epsilon__dot^2)*sin(epsilon) - 2*r__a*(r__a*alpha__hatdot*theta__dot - s__2*alpha__hatdot*theta__dot)*sin(alpha__hat);
equation_2=F_equation_2(alpha__hat,theta,epsilon,alpha__hatdot,theta__dot,epsilon__dot,s__1,s__2);

F_equationextra_2=@(alpha__hat,theta,epsilon,s__1,s__2) -v*cos(theta)*(theta__dot)^2*F_xddotcoeff_2(alpha__hat,theta,epsilon,s__1,s__2) -v*sin(theta)*(theta__dot)^2*F_yddotcoeff_2(alpha__hat,theta,epsilon,s__1,s__2);
equationextra_2=F_equationextra_2(alpha__hat,theta,epsilon,s__1,s__2);

F_equation_2_final=@(alpha__hat,theta,epsilon,alpha__hatdot,theta__dot,epsilon__dot,s__1,s__2) F_equation_2(alpha__hat,theta,epsilon,alpha__hatdot,theta__dot,epsilon__dot,s__1,s__2)+F_equationextra_2(alpha__hat,theta,epsilon,s__1,s__2);
equation_2_final=F_equation_2_final(alpha__hat,theta,epsilon,alpha__hatdot,theta__dot,epsilon__dot,s__1,s__2);

F_Q_2=@(theta,epsilon,epsilon__dot) I_w*cos(epsilon)*epsilon__dot;%I_w*w*(sin(theta)*cos(epsilon)*theta + cos(theta)*sin(epsilon)*epsilon__dot);
Q_2=F_Q_2(theta,epsilon,epsilon__dot);
end