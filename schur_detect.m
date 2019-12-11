clear all
clc
A = [0 1 0;
     1 1 0;
     -1 0 0];
 B = [0;1;0];
 C=[1 0 0;
     0 0 1];
 P = sdpvar(3,3);
 W = sdpvar(3,2);
 gamma = sdpvar(1);
 M11 = -P;
 M12 = A'*P + C'*W';
 M21 = P*A+W*C;
 M22 = -P;
 M = [M11 M12;
     M21 M22];
 N =[-P P*A;
     A'*P -P-gamma*C'*C];
 F=[];
 F=[F;M<zeros(size(M))];
 F=[F;N,zeros(size(N))];
 optimize(F);
 P=value(P);
 eig(P)
 %The system is not Schur detectable