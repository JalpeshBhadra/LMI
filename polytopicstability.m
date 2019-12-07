%Polytopic Stability
clear all
clc

%Given Variables
tau_z = 1;
Jxyz = 0.7501;
g = 9.8;
V = 1401;
b7 = -.001;
wx = 10;


%From Table 11.2 Duang
a1t = [1.593	1.485	1.269	1.13	0.896	0.559	0.398];
ap1t= [0.285	0.192	0.147	0.118	0.069	0.055	0.043];
a2t = [260.559	266.415	196.737	137.385	129.201	66.338	51.003];
a3t = [185.488	182.532	176.932	160.894	138.591	78.404	53.84];
a4t = [1.506	1.295	1.169	1.13	1.061	0.599	0.421];
a5t = [0.298	0.243	0.217	0.191	0.165	0.105	0.078];


%Variables 
P = sdpvar(3,3);

tol = 0.00001;

for i=1:1:7
A=[-a4t(i) 1 -a5t(i);((-ap1t(i)*a4t(i))-a2t(i)) (ap1t(i)-a1t(i)) ((-ap1t(i)*a5t(i))-a3t(i));0 0 -(1/tau_z)];



end
options=sdpsettings('solver','sedumi');
F=[];
F=[F; (A'*P+P*A)<-tol*eye(3)];
optimize(F);
