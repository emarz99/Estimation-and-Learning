function [A,B,C,D] = quad2grey(Xu,Xq,Mu,Mq,Xd,Md,Ts)
 
A=[Xu, Xq, -9.81; 
    Mu, Mq, 0;
    0, 1, 0];


B=[Xd; Md; 0];

C=[ 0, 1, 0;
    Xu, Xq, 0]; 


D=[ 0; Xd];
end