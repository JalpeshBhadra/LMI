clear all
clc
A= [-0.0558 -0.9968 0.0802 0.0415;
    0.598 -0.1150 -0.0318 0;
    -3.05 0.388 -0.465 0;
    0 0.08 1 0];
B=[0.0729 0.0001;
    -4.75 1.23;
    1.53 10.63;
    0 0];
C=[0 1 0 0;
    0 0 0 1];
%% 
gamma = sdpvar(1);
W = sdpvar(2,4);
P = sdpvar(4,4);
F1=[];
F1=[F1;A'*P+P*A+W'*C+C'*W < zeros(4,4)];
F1=[F1;A'*P+P*A-gamma*C'*C < zeros(4,4)];
F1=[F1;C*(A'*P+P*A)*C' < zeros(2,2)];
options = sdpsettings('solver','sedumi');
optimize(F1);
gamma = value(gamma)
P = value(P)
eig(P)
%System is not Hurwitz Detectable
% System is not 