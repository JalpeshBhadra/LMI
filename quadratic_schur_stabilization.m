%Quadratic schur stabilization
A=[0.8147*10    0.9134*10    0.2785*10;
   0.9058*5    0.6324*20    0.5469*10;
   0.1270*10    0.0975*10    0.9575*15]
Ai = [0.9649    0.9572    0.1419;
     0.1576    0.4854    0.4218;
     0.9706    0.8003    0.9157]
 B = [0.7922*10;
    0.9595*20;
    0.6557*15]
Bi = [0.0357;
    0.8491;
    0.9340]
X= sdpvar(3,3)
Z= sdpvar(1,3)
tol = 0.00001;
M = [X A*X+B*Z ; (A*X+B*Z)' X]
N = [zeros(3) Ai*X+Bi*Z ; (Ai*X+Bi*Z)' zeros(3)]

F=[F;X > zeros(3,3)]
F=[F;M+N > tol*eye(size(M))]

XX = value(X)
ZZ = value(Z)
% Extracting controller
K = ZZ*(inv(XX))
