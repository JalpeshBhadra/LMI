%H2 optimal state feedback controller
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
Y = sdpvar(3,3);
Z = sdpvar(1,3);
gamma = sdpvar(1);
F = [];

tol = 0.00001;
W = eye(3)
C1 = eye(3)
for i=1:1:7
A=[-a4t(i) 1 -a5t(i);((-ap1t(i)*a4t(i))-a2t(i)) (ap1t(i)-a1t(i)) ((-ap1t(i)*a5t(i))-a3t(i));0 0 -(1/tau_z)];
B1 = (wx/57.3)*[-1 0;-ap1t(i) Jxyz; 0 0];
B2 = [0;0;(1/tau_z)]; 
C = (1/(57.3*g))*[(57.3*g) 0 0;V*a4t(i) 0 V*a5t(i)];
D2 = [0;0]; %Add 0 as 2nd row (Correct in D)
D1 = (1/(57.3*g))*[0 0; V*b7 0];

M11 = A*Y + B2*Z + Y*A' + Z'*B2' 
M12 = B1
M21 = B1'
M22 = -eye(2)
M = [M11 M12; M21 M22]

N11 = Y 
N12 = (C1*Y + zeros(3))'
N21 = (C1*Y + zeros(3))
N22 = W

N = [N11 N12; N21 N22]

trace(W) < -eye(3)*gamma;

end
options=sdpsettings('solver','sedumi');
F=[F;M<tol*eye(size(M))];
F= [F;N>tol*eye(size(N))];
optimize(F,gamma,options);

value(gamma)
%Controller
K = value(Z)*(value(Y'))
