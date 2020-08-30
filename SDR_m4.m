 %Implementing Semi-definite programming method for m = 4 %Compare the result of different constraints
clear;
format long
n=5;m=4;In=eye(n);Im=eye(m);     
%Generate a random system G(s)
A=-abs(diag(randn(n,1)));
B=randn(n,1);
C=randn(1,n);
G = ss(A,B,C,0);

%Implement Semi-definite program for m = 4
[g1,L1,LL1,Gm2_1,pp1]=SDR_forth(A,B,C);
[g2,L2,LL2,Gm2_2,pp2]=SDR_forth_sedumi(A,B,C);
lamda1=min(eig(L1));
lamda2=min(eig(L2));
%Compute the interpolating points
ss1=p2s([pp1(1),pp1(n+1),pp1(2*n+1),pp1(3*n+1)]); 
ss2=p2s([pp2(1),pp2(n+1),pp2(2*n+1),pp2(3*n+1)]); 
%IRKA
[Gm_min1,iter_no1,sigma1] = IRKA(n,m,G,ss1); 
[Gm_min2,iter_no2,sigma2] = IRKA(n,m,G,ss2); 
Gm1=Gm_generate(n,m,G,ss1);
Gm2=Gm_generate(n,m,G,ss2);
% bode(G,Gm1,Gm2);
% set(findobj(gcf,'type','line'),'linewidth',1);

function [g,L,LL,Gm2,pp]=SDR_forth(A,B,C) 
n=size(A,1);m=4;I=eye(n);O=0*I;
Ai=inv(A);Ai2=Ai*Ai;Ai3=Ai2*Ai;Ai4=Ai3*Ai; BB=B*B';
T6=[Ai*BB*Ai';
    O;
Ai*BB*Ai3'+Ai3*BB*Ai'-Ai2*BB*Ai2';
O;
O]; 
T7=[O,Ai2*BB*Ai2',O,Ai4*BB*Ai2'+Ai2*BB*Ai4'-Ai3*BB*Ai3',O;
O,O,Ai3*BB*Ai3',O,O;
O,O,O,Ai4*BB*Ai4',O;
O,O,O,O,O;
O,O,O,O,O];
T4=[Ai;
    -Ai2;
    Ai3;
    -Ai4;
    O];
cvx_begin sdp; 
cvx_precision best
cvx_solver SDPT3
variable g(1);
variable S1(n,n);
variable S2(n,n);
variable S3(n,n);
variable S4(n,n);
variable S5(n,n);
variable G12(n,n);
variable G13(n,n);
variable G14(n,n);
variable G23(n,n);
variable G24(n,n);
variable G34(n,n);
variable P(n,n);
variable Q(n,n);
minimize(g)
subject to 
S=[S1;S2;S3+Q;S4;S5];
G=[O,G12-Q,G13,G14-P,Q; 
   (G12-Q)',O,G23,G24,O;
G13',G23',O,G34,P;
(G14-P)',G24',G34',O,O;
(Q)',O,(P)',O,O]; 

S1+S1'>=0;
S2+S2'>=0;
S3+S3'>=0;
S4+S4'>=0; 
S5+S5'>=0;
G12+G12'>=0;
G13+G13'>=0;
G14+G14'>=0;
G23+G23'>=0;
G24+G24'>=0;
G34+G34'>=0;
P+P'>=0;

[g,C*(T6+S)'; (T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]>=0;

cvx_end

S1=full(S1);S2=full(S2);S3=full(S3);S4=full(S4);S5=full(S5);
P=full(P);Q=full(Q);
G12=full(G12);G13=full(G13);G14=full(G14);
G23=full(G23);G24=full(G24);G34=full(G34);
S=[S1;S2;S3+Q;S4;S5];
G=[O,G12-Q,G13,G14-P,Q; 
   (G12-Q)',O,G23,G24,O;
G13',G23',O,G34,P;
(G14-P)',G24',G34',O,O;
(Q)',O,(P)',O,O]; 
L=full([g,C*(T6+S)'; (T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]);
LL=L;
L=L(:,1:1+m*n);
L=L(1:1+m*n,:);
%Construct vector pp
[t,tt]=eig(L);[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
T6=[Ai*BB*Ai'; O;
Ai*BB*Ai3'+Ai3*BB*Ai'-Ai2*BB*Ai2';
O]; 
T7=[O,Ai2*BB*Ai2',O,Ai4*BB*Ai2'+Ai2*BB*Ai4'-Ai3*BB*Ai3';
O,O,Ai3*BB*Ai3',O; O,O,O,Ai4*BB*Ai4'; O,O,O,O];
T4=[Ai; -Ai2;Ai3; -Ai4];
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((Gm2-g)/Gm2);

x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);x3=x(2*n+1:3*n);x4=x(3*n+1:4*n); 
z=(T4'*x-C'); 
pp1=x1./z; 
pp2=x2./z; 
pp3=x3./z; 
pp4=x4./z; 
pp=[pp1;pp2;pp3;pp4];
end

function [g,L,LL,Gm2,pp]=SDR_forth_sedumi(A,B,C) 
n=size(A,1);m=4;I=eye(n);O=0*I;
Ai=inv(A);Ai2=Ai*Ai;Ai3=Ai2*Ai;Ai4=Ai3*Ai; BB=B*B';
T6=[Ai*BB*Ai';
    O;
Ai*BB*Ai3'+Ai3*BB*Ai'-Ai2*BB*Ai2';
O;
O]; 
T7=[O,Ai2*BB*Ai2',O,Ai4*BB*Ai2'+Ai2*BB*Ai4'-Ai3*BB*Ai3',O;
O,O,Ai3*BB*Ai3',O,O;
O,O,O,Ai4*BB*Ai4',O;
O,O,O,O,O;
O,O,O,O,O];
T4=[Ai;
    -Ai2;
    Ai3;
    -Ai4;
    O];
cvx_begin sdp; 
cvx_precision best
cvx_solver sedumi
variable g(1);
variable S1(n,n);
variable S2(n,n);
variable S3(n,n);
variable S4(n,n);
variable S5(n,n);
variable G12(n,n);
variable G13(n,n);
variable G14(n,n);
variable G23(n,n);
variable G24(n,n);
variable G34(n,n);
variable P(n,n);
variable Q(n,n);
minimize(g)
subject to 
S=[S1;S2;S3+Q;S4;S5];
G=[O,G12-Q,G13,G14-P,Q; 
   (G12-Q)',O,G23,G24,O;
G13',G23',O,G34,P;
(G14-P)',G24',G34',O,O;
(Q)',O,(P)',O,O]; 

S1+S1'>=0;
S2+S2'>=0;
S3+S3'>=0;
S4+S4'>=0; 
S5+S5'>=0;
G12+G12'>=0;
G13+G13'>=0;
G14+G14'>=0;
G23+G23'>=0;
G24+G24'>=0;
G34+G34'>=0;
P+P'>=0;

[g,C*(T6+S)'; (T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]>=0;

cvx_end

S1=full(S1);S2=full(S2);S3=full(S3);S4=full(S4);S5=full(S5);
P=full(P);Q=full(Q);
G12=full(G12);G13=full(G13);G14=full(G14);
G23=full(G23);G24=full(G24);G34=full(G34);
S=[S1;S2;S3+Q;S4;S5];
G=[O,G12-Q,G13,G14-P,Q; 
   (G12-Q)',O,G23,G24,O;
G13',G23',O,G34,P;
(G14-P)',G24',G34',O,O;
(Q)',O,(P)',O,O]; 
L=full([g,C*(T6+S)'; (T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]);
LL=L;
L=L(:,1:1+m*n);
L=L(1:1+m*n,:);
%Construct vector pp
[t,tt]=eig(L);[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
T6=[Ai*BB*Ai'; O;
Ai*BB*Ai3'+Ai3*BB*Ai'-Ai2*BB*Ai2';
O]; 
T7=[O,Ai2*BB*Ai2',O,Ai4*BB*Ai2'+Ai2*BB*Ai4'-Ai3*BB*Ai3';
O,O,Ai3*BB*Ai3',O; O,O,O,Ai4*BB*Ai4'; O,O,O,O];
T4=[Ai; -Ai2;Ai3; -Ai4];
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((Gm2-g)/Gm2);

x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);x3=x(2*n+1:3*n);x4=x(3*n+1:4*n); 
z=(T4'*x-C'); 
pp1=x1./z; 
pp2=x2./z; 
pp3=x3./z; 
pp4=x4./z; 
pp=[pp1;pp2;pp3;pp4];
end

function [ss]=p2s(p)
syms s;
ss=solve(s^4-p(1)*s^3+p(2)*s^2-p(3)*s+p(4)==0,s);
ss=double(ss);
end
