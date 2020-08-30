% Demonstration of Semi-Definite Relaxation Method for m = 2 %n = 3, m = 2
clear;
format long;

% 
%Generate a random system G(s) 
n=35;m=2;In=eye(n);Im=eye(m);
A=-diag(abs(randn(n,1))); 
B=randn(n,1);
C=randn(1,n);

G = ss(A,B,C,0);


% examples in previous literature

% numerator=[-1.3369 -4.8341 -47.5819 -42.7285];
% denominator=[1 17.0728 84.9908 122.4400 59.9309];
% sys = tf(numerator,denominator); 
% [A,B,C,D] = tf2ss(numerator,denominator); 
% G=ss(A,B,C,D);


n=size(A,1);m=2;
In=eye(n);Im=eye(m);

%Implement the SDR method
[g1,L1,S1_1,S2_1,G12_1,Gm2_1,pp1]=SDR_second(A,B,C); 
% [g1,L1,S1_1,S2_1,G12_1,Gm2_1,pp1]=SDR_second_modified(A,B,C); 
[g2,L2,S1_2,S2_2,G12_2,Gm2_2,pp2]=SDR_second_sedumi(A,B,C);
%Results
lamda1 = min(eig(L1));
lamda2 = min(eig(L2));
%Calculate interpolating points
syms s;
ss1=solve(s^2-pp1(1)*s+pp1(n+1)==0,s);
ss2=solve(s^2-pp2(1)*s+pp2(n+1)==0,s);
ss1=double(ss1);
ss2=double(ss2);
%Implement the IRKA method
[Gm_min1,iter_no1,sigma1] = IRKA(n,m,G,ss1);
[Gm_min2,iter_no2,sigma2] = IRKA(n,m,G,ss2);
%bode plot
% bode(G,Gm_min1,Gm_min2);
% set(findobj(gcf,'type','line'),'linewidth',1);


function [g,L,S1,S2,G12,Gm2,pp]=SDR_second(A,B,C) 
n=size(A,1);m=2;I=eye(n);O=0*I;
Ai=inv(A);Ai2=Ai*Ai;
T6=[Ai*B*B'*Ai';O]; 
T7=[O,Ai2*B*B'*Ai2';O,O]; 
T4=[Ai;-Ai2];

cvx_begin sdp; 
cvx_precision best 
variable g(1); 
variable S1(n,n); 
variable S2(n,n); 
variable G12(n,n); 
minimize(g) 
subject to 
S=[S1;S2]; 
G=[O,G12;G12',O]; 
S1+S1'>=O;
S2+S2'>=O; 
G12+G12'>=0; 
[g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]>=0; 
cvx_end

S1=full(S1);S2=full(S2);G12=full(G12);  
S=[S1;S2];
G=[O,G12;G12',O]; 
L=full([g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]);
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((g-Gm2)/Gm2);
x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);
z=(T4'*x-C'); 
pp1=x1./z;
pp2=x2./z;
pp=[pp1;pp2];
end

function [g,L,S1,S2,G12,Gm2,pp]=SDR_second_modified(A,B,C) 
n=size(A,1);m=2;I=eye(n);O=0*I;
Ai=inv(A);Ai2=Ai*Ai;
T6=[Ai*B*B'*Ai';O]; 
T7=[O,Ai2*B*B'*Ai2';O,O]; 
T4=[Ai;-Ai2];

cvx_begin sdp; 
cvx_precision best 
variable g(1); 
variable S1(n,n); 
variable S2(n,n); 
variable G12(n,n); 
minimize(g) 
subject to 
S=[S1;S2]; 
G=[O,G12;G12',O]; 
S1+S1'>=O;
S2+S2'>=O; 
G12+G12'==0; 
[g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]>=0; 
cvx_end

S1=full(S1);S2=full(S2);G12=full(G12);  
S=[S1;S2];
G=[O,G12;G12',O]; 
L=full([g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]);
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((g-Gm2)/Gm2);
x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);
z=(T4'*x-C'); 
pp1=x1./z;
pp2=x2./z;
pp=[pp1;pp2];
end

function [g,L,S1,S2,G12,Gm2,pp]=SDR_second_sedumi(A,B,C) 
n=size(A,1);m=2;I=eye(n);O=0*I;
Ai=inv(A);Ai2=Ai*Ai;
T6=[Ai*B*B'*Ai';O]; 
T7=[O,Ai2*B*B'*Ai2';O,O]; 
T4=[Ai;-Ai2];

cvx_begin sdp; 
cvx_precision best 
cvx_solver sedumi
variable g(1); 
variable S1(n,n); 
variable S2(n,n); 
variable G12(n,n); 
minimize(g) 
subject to 
S=[S1;S2]; 
G=[O,G12;G12',O]; 
S1+S1'>=O;
S2+S2'>=O; 
G12+G12'>=0; 
[g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]>=0; 
cvx_end

S1=full(S1);S2=full(S2);G12=full(G12);  
S=[S1;S2];
G=[O,G12;G12',O]; 
L=full([g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]);
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((g-Gm2)/Gm2);
x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);
z=(T4'*x-C'); 
pp1=x1./z;
pp2=x2./z;
pp=[pp1;pp2];
end

 