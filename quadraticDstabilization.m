

A = [-1.3410 0.9933 0 -0.1689 -0.2518
    43.2230 -0.8693 0 -17.2510 -1.5766
    1.3410 0.0067 0 0.1689 0.2518
    0 0 0 -20.0000 0
    0 0 0 0 -20.0000];
B = [0 0
    0 0
    0 0
    20 0
    0 20];
Ai =[0.6787    0.1712    0.0971    0.0344    0.1869;
    0.7577    0.7060    0.8235    0.4387    0.4898;
    0.7431    0.0318    0.6948    0.3816    0.4456;
    0.3922    0.2769    0.3171    0.7655    0.6463;
    0.6555    0.0462    0.9502    0.7952    0.7094]
Bi = [ 0.7547    0.1190;
    0.2760    0.4984;
    0.6797    0.9597;
    0.6551    0.3404;
    0.1626    0.5853]
%matrices that will be used to demonstrate the controller
[~,n]=size(A);[c1,c2]=size(B);
Z=sdpvar(c2,c1);P3=sdpvar(n,n);%initialize the Z and P 
%sdpvar matrices for this problem.
options = sdpsettings('solver','sedumi');
F1=[P3>=1e-5*eye(n)];
ts=6;p_os=0.1;tr=1;r=(1.8/tr);alpha=(4.6/ts);c=(log(p_os)/pi);
M = [(-r*P3) (A*P3+B*Z); (A*P3+B*Z)' (-r*P3)]
N = [ zeros(5,5) (Ai*P3+Bi*Z);(Ai*P3+Bi*Z)' zeros(5,5)]
F1=[F1,(M+N)<=1e-5*eye(2*n)];
F1=[F1,A*P3+B*Z+(A*P3+B*Z)'+2*alpha*P3+Ai*P3+Bi*Z+(Ai*P3+Bi*Z)'<=1e-5*eye(n)];
R = [A*P3+B*Z+(A*P3+B*Z)' c*(A*P3+B*Z-(A*P3+B*Z)'); c*((A*P3+B*Z)'-(A*P3+B*Z)) A*P3+B*Z+(A*P3+B*Z)']
W = [Ai*P3+Bi*Z+(Ai*P3+Bi*Z)' c*(Ai*P3+Bi*Z-(Ai*P3+Bi*Z)'); c*((Ai*P3+Bi*Z)'-(Ai*P3+Bi*Z)) Ai*P3+Bi*Z+(Ai*P3+Bi*Z)']
F1=[F1,R+W<=1e-5*eye(2*n)];
optimize(F1,[],options);%execute the given problem
Z=value(Z);P3=value(P3);
K=Z*inv(P3)%display the resulting controller K matrix
ApBK_eig=eig(A+B*K)%display the resulting eigenvalues of A+BK
%Note that all the real part of eigenvalues of A+BK are negative, hence
%stable.