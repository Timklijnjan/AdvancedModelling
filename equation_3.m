function [alphaddotcoeff_3,thetaddotcoeff_final_3,epsddotcoeff_3,Q_3,equation_3_final]=equation_3(alpha__hat,theta,epsilon,alpha__hatdot,theta__dot,epsilon__dot,v,g)
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
alphaddotcoeff_3=I__xy;

F_thetaddotcoeff_3=@(alpha__hat,epsilon,s__1,s__2) -(2*r__a*s__1 - 2*s__1*s__2)*cos(epsilon)*M__a*cos(alpha__hat)^2 + ((-r__a^2 + 2*r__a*s__2 + s__1^2 - s__2^2)*sin(alpha__hat) + r__a*s__1)*cos(epsilon)*M__a*cos(alpha__hat) - M__a*((-r__a^2 + r__a*s__2)*sin(alpha__hat) - r__a*s__1 + s__1*s__2)*cos(epsilon) + I__xz;
thetaddotcoeff_3=F_thetaddotcoeff_3(alpha__hat,epsilon,s__1,s__2);

F_epsddotcoeff_3=@(alpha__hat,s__1,s__2) -(-r__a^2 + 2*r__a*s__2 + s__1^2 - s__2^2)*M__a*cos(alpha__hat)^2 + ((-2*r__a*s__1 + 2*s__1*s__2)*sin(alpha__hat) + 2*r__a*(s__2 - r__a))*M__a*cos(alpha__hat) + 2*r__a*s__1*M__a*sin(alpha__hat) + (r__a^2 + s__1^2)*M__a + I__xx;
epsddotcoeff_3=F_epsddotcoeff_3(alpha__hat,s__1,s__2);

F_xddotcoeff_3=@(alpha__hat,theta,epsilon,s__1,s__2) (-s__2 + r__a)*sin(theta)*cos(epsilon)*M__a*cos(alpha__hat) - M__a*(sin(theta)*s__1*sin(alpha__hat) + sin(theta)*r__a)*cos(epsilon);
xddotcoeff_3=F_xddotcoeff_3(alpha__hat,theta,epsilon,s__1,s__2);

F_yddotcoeff_3=@(alpha__hat,theta,epsilon,s__1,s__2) -(-s__2 + r__a)*cos(theta)*cos(epsilon)*M__a*cos(alpha__hat) - M__a*(-cos(theta)*s__1*sin(alpha__hat) - cos(theta)*r__a)*cos(epsilon);
yddotcoeff_3=F_yddotcoeff_3(alpha__hat,theta,epsilon,s__1,s__2);

F_thetaddotextra_x_3=@(alpha__hat,theta,epsilon,s__1,s__2)-v*sin(theta)*F_xddotcoeff_3(alpha__hat,theta,epsilon,s__1,s__2);
thetaddotextra_x_3=F_thetaddotextra_x_3(alpha__hat,theta,epsilon,s__1,s__2);

F_thetaddotextra_y_3=@(alpha__hat,theta,epsilon,s__1,s__2) v*cos(theta)*F_yddotcoeff_3(alpha__hat,theta,epsilon,s__1,s__2);
thetaddotextra_y_3=F_thetaddotextra_y_3(alpha__hat,theta,epsilon,s__1,s__2);

F_thetaddotcoeff_final_3=@(alpha__hat,theta,epsilon,s__1,s__2) F_thetaddotcoeff_3(alpha__hat,epsilon,s__1,s__2) %mistake x=vcos(theta): +F_thetaddotextra_x_3(alpha__hat,theta,epsilon,s__1,s__2)+F_thetaddotextra_y_3(alpha__hat,theta,epsilon,s__1,s__2);
thetaddotcoeff_final_3=F_thetaddotcoeff_final_3(alpha__hat,theta,epsilon,s__1,s__2);

F_equation_3=@(alpha__hat,theta,epsilon,s__1,s__2,theta__dot,alpha__hatdot,epsilon__dot) -((theta__dot^2*(r__a + s__1 - s__2)*(r__a - s__1 - s__2)*sin(epsilon) + 2*r__a^2*alpha__hatdot*theta__dot - 4*r__a*alpha__hatdot*theta__dot*s__2 - 2*alpha__hatdot*theta__dot*s__1^2 + 2*alpha__hatdot*theta__dot*s__2^2)*cos(epsilon) + 4*r__a*s__1*alpha__hatdot*epsilon__dot - 4*s__1*alpha__hatdot*epsilon__dot*s__2)*M__a*cos(alpha__hat)^2 + (((2*theta__dot^2*s__1*(-s__2 + r__a)*sin(epsilon) + 4*r__a*s__1*alpha__hatdot*theta__dot - 4*s__1*alpha__hatdot*theta__dot*s__2)*sin(alpha__hat) + 2*theta__dot^2*r__a*(-s__2 + r__a)*sin(epsilon) + 2*r__a*(r__a*alpha__hatdot*theta__dot - s__2*alpha__hatdot*theta__dot))*cos(epsilon) + (-2*r__a^2*alpha__hatdot*epsilon__dot + 4*r__a*s__2*alpha__hatdot*epsilon__dot + 2*s__1^2*alpha__hatdot*epsilon__dot - 2*s__2^2*alpha__hatdot*epsilon__dot)*sin(alpha__hat) + g*(-s__2 + r__a)*sin(epsilon) + 2*r__a*s__1*alpha__hatdot*epsilon__dot)*M__a*cos(alpha__hat) - M__a*((2*theta__dot^2*sin(epsilon)*r__a*s__1 + 2*r__a*s__1*alpha__hatdot*theta__dot)*sin(alpha__hat) + theta__dot^2*(r__a^2 + s__1^2)*sin(epsilon) + 2*alpha__hatdot*theta__dot*s__1^2)*cos(epsilon) - (g*sin(epsilon)*s__1 - 2*r__a*(r__a*alpha__hatdot*epsilon__dot - s__2*alpha__hatdot*epsilon__dot))*M__a*sin(alpha__hat) - M__a*sin(epsilon)*g*r__a + (2*r__a*s__1*alpha__hatdot*epsilon__dot - 2*s__1*s__2*alpha__hatdot*epsilon__dot)*M__a;
equation_3=F_equation_3(alpha__hat,theta,epsilon,s__1,s__2,theta__dot,alpha__hatdot,epsilon__dot);

F_equationextra_3=@(alpha__hat,theta,epsilon,s__1,s__2) -v*sin(theta)*theta__dot*F_xddotcoeff_3(alpha__hat,theta,epsilon,s__1,s__2) +v*cos(theta)*theta__dot*F_yddotcoeff_3(alpha__hat,theta,epsilon,s__1,s__2);%mistake x=vcos(theta): -v*cos(theta)*(theta__dot)^2*F_xddotcoeff_3(alpha__hat,theta,epsilon,s__1,s__2) -v*sin(theta)*(theta__dot)^2*F_yddotcoeff_3(alpha__hat,theta,epsilon,s__1,s__2);
equationextra_3=F_equationextra_3(alpha__hat,theta,epsilon,s__1,s__2);

F_equation_3_final=@(alpha__hat,theta,epsilon,s__1,s__2,theta__dot,alpha__hatdot,epsilon__dot) F_equation_3(alpha__hat,theta,epsilon,s__1,s__2,theta__dot,alpha__hatdot,epsilon__dot)+F_equationextra_3(alpha__hat,theta,epsilon,s__1,s__2);
equation_3_final=F_equation_3_final(alpha__hat,theta,epsilon,s__1,s__2,theta__dot,alpha__hatdot,epsilon__dot);

F_Q_3= @(theta,epsilon,epsilon__dot,theta__dot) I_w*w*(cos(theta)*cos(epsilon)*theta__dot - sin(theta)*sin(epsilon)*epsilon__dot);%@(epsilon) I_w*w*epsilon__dot*cos(epsilon);
Q_3=F_Q_3(theta,epsilon,epsilon__dot,theta__dot);
end