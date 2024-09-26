function [A,B,C,D] = quadgrey(Xu,Xq,Mu,Mq,Xd,Md,Ts)
 
A=[Xu, Xq, -9.81; 
    Mu, Mq, 0;
    0, 1, 0];


B=[Xd; Md; 0];

C=[1, 0, 0;
    0, 1, 0;
    0, 0, 1;
    Xu, Xq, 0]; 


D=[0;0; 0; Xd];
end
