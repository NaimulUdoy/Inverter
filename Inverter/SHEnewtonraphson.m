clc;
clear all;
A=input('M='); 
H=input('H=');
M=0:0.01:A; 

for ii = 1:length(M)

    alpha1=20*pi/180;
    alpha2=30*pi/180;
    alpha3=50*pi/180;
    alpha4=70*pi/180;
    alpha5=80*pi/180;

 for i=1:H
       
        T = [M(ii)*pi/4 0 0 0 0]';

        % F matrix: Cosine terms for harmonics
        F = [cos(alpha1)-cos(alpha2)+cos(alpha3)-cos(alpha4)+cos(alpha5);
            cos(3*alpha1)-cos(3*alpha2)+cos(3*alpha3)-cos(3*alpha4)+cos(3*alpha5);
            cos(5*alpha1)-cos(5*alpha2)+cos(5*alpha3)-cos(5*alpha4)+cos(5*alpha5);
            cos(7*alpha1)-cos(7*alpha2)+cos(7*alpha3)-cos(7*alpha4)+cos(7*alpha5);
            cos(9*alpha1)-cos(9*alpha2)+cos(9*alpha3)-cos(9*alpha4)+cos(9*alpha5)];

      dF =[-sin(alpha1) +sin(alpha2) -sin(alpha3) +sin(alpha4) -sin(alpha5);
           -3*sin(3*alpha1) +3*sin(3*alpha2) -3*sin(3*alpha3)  +3*sin(3*alpha4) -3*sin(3*alpha5);
           -5*sin(5*alpha1) +5*sin(5*alpha2) -5*sin(5*alpha3) +5*sin(5*alpha4) -5*sin(5*alpha5);
           -7*sin(7*alpha1) +7*sin(7*alpha2) -7*sin(7*alpha3) +7*sin(7*alpha4) -7*sin(7*alpha5);
           -9*sin(9*alpha1) +9*sin(9*alpha2) -9*sin(9*alpha3) +9*sin(9*alpha4) -9*sin(9*alpha5)];

        dalpha=(inv(dF))*(T-F);

       
        i;
        alpha=[alpha1;alpha2;alpha3;alpha4;alpha5]*180/pi

        % Update alpha values with the calculated delta
        F;
        dalpha*180/pi;
        alpha1=alpha1+dalpha(1);
        alpha2=alpha2+dalpha(2);
        alpha3=alpha3+dalpha(3);
        alpha4=alpha4+dalpha(4);
        alpha5=alpha5+dalpha(5);

        % Check for convergence
        if norm(dalpha)<1e-15
        break;
        end 
end
end

