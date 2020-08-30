clear;
format long
n=4;m=3;In=eye(n);Im=eye(m);  
%Generate a random system G(s)
A=-abs(diag(randn(n,1)));
B=randn(n,1);
C=randn(1,n);
G = ss(A,B,C,0);

%Implement the SDR method
[g,L,S1,S2,S3,G3,G12,G13,G23,Gm2_or,pp]=SDR_third(A,B,C) ;
[g_s,L_s,S1_s,S2_s,S3_s,G3_s,G12_s,G13_s,G23_s,Gm2_or_s,pp_s]=SDR_third_sedumi(A,B,C) ;

[g1,L1,LL1,Gm2_al,pp1]=SDR_third_alternative(A,B,C);
[g1_s,L1_s,LL1_s,Gm2_al_s,pp1_s]=SDR_third_alternative_sedumi(A,B,C);
ss1=p2s([pp1(n),pp1(2*n),pp1(3*n)]); 
ss=p2s([pp(n),pp(2*n),pp(3*n)]); 
ss1_s=p2s([pp1_s(n),pp1_s(2*n),pp1_s(3*n)]); 
ss_s=p2s([pp_s(n),pp_s(2*n),pp_s(3*n)]); 
%Run IRKA
[Gm_min1,iter_no1,sigma1] = IRKA(n,m,G,ss1); 
[Gm_min,iter_no,sigma] = IRKA(n,m,G,ss);
[Gm_min1_s,iter_no1_s,sigma1_s] = IRKA(n,m,G,ss1_s); 
[Gm_min_s,iter_no_s,sigma_s] = IRKA(n,m,G,ss_s);
lamda1=min(eig(L1));
lamda=min(eig(L));
lamda1_s=min(eig(L1_s));
lamda_s=min(eig(L_s));
Gm1=Gm_generate(n,m,G,ss1);Gm=Gm_generate(n,m,G,ss);
Gm1_s=Gm_generate(n,m,G,ss1_s);Gm_s=Gm_generate(n,m,G,ss_s);

% bode(G,Gm,Gm_s,Gm1,Gm1_s);
% set(findobj(gcf,'type','line'),'linewidth',1);

function [g,L,S1,S2,S3,G3,G12,G13,G23,Gm2,pp]=SDR_third(A,B,C) 
n=size(A,1);m=3;
I=eye(n);O=0*I;Ai=inv(A);Ai2=Ai*Ai;Ai3=Ai2*Ai;BB=B*B'; 
T6=[Ai*BB*Ai';O;Ai*BB*Ai3'+Ai3*BB*Ai'-Ai2*BB*Ai2']; 
T7=[O,Ai2*BB*Ai2',O;O,O,Ai3*BB*Ai3';O,O,O]; 
T4=[Ai;-Ai2;Ai3];

cvx_begin sdp; 
cvx_precision best
variable g(1); 
variable S1(n,n); 
variable S2(n,n); 
variable S3(n,n); 
variable G3(n,n); 
variable G12(n,n); 
variable G13(n,n); 
variable G23(n,n); 
minimize(g) 
subject to 
S=[S1;S2;S3-G3]; 
G=[O,G12+G3,G13;G12'+G3',O,G23;G13',G23',O]; 
S1+S1'>=O;
S2+S2'>=O; 
S3+S3'>=O;
G3+G3'>=O;
G12+G12'>=0;
G13+G13'>=0;
G23+G23'>=0;
[g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]>=0; 
cvx_end

S1=full(S1);S2=full(S2);S3=full(S3);G3=full(G3);
G12=full(G12);G13=full(G13);G23=full(G23);
S=[S1;S2;S3-G3]; 
G=[O,G12+G3,G13;G12'+G3',O,G23;G13',G23',O]; 
L=full([g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]);
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((Gm2-g)/Gm2); 
x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);x3=x(2*n+1:3*n);
z=(T4'*x-C'); 
pp1=x1./z;
pp2=x2./z;
pp3=x3./z;
pp=[pp1;pp2;pp3];
end

function [g,L,S1,S2,S3,G3,G12,G13,G23,Gm2,pp]=SDR_third_sedumi(A,B,C) 
n=size(A,1);m=3;
I=eye(n);O=0*I;Ai=inv(A);Ai2=Ai*Ai;Ai3=Ai2*Ai;BB=B*B'; 
T6=[Ai*BB*Ai';O;Ai*BB*Ai3'+Ai3*BB*Ai'-Ai2*BB*Ai2']; 
T7=[O,Ai2*BB*Ai2',O;O,O,Ai3*BB*Ai3';O,O,O]; 
T4=[Ai;-Ai2;Ai3];

cvx_begin sdp; 
cvx_precision best
cvx_solver sedumi
variable g(1); 
variable S1(n,n); 
variable S2(n,n); 
variable S3(n,n); 
variable G3(n,n); 
variable G12(n,n); 
variable G13(n,n); 
variable G23(n,n); 
minimize(g) 
subject to 
S=[S1;S2;S3-G3]; 
G=[O,G12+G3,G13;G12'+G3',O,G23;G13',G23',O]; 
S1+S1'>=O;
S2+S2'>=O; 
S3+S3'>=O;
G3+G3'>=O;
G12+G12'>=0;
G13+G13'>=0;
G23+G23'>=0;
[g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]>=0; 
cvx_end

S1=full(S1);S2=full(S2);S3=full(S3);G3=full(G3);
G12=full(G12);G13=full(G13);G23=full(G23);
S=[S1;S2;S3-G3]; 
G=[O,G12+G3,G13;G12'+G3',O,G23;G13',G23',O]; 
L=full([g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G]);
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((Gm2-g)/Gm2); 
x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);x3=x(2*n+1:3*n);
z=(T4'*x-C'); 
pp1=x1./z;
pp2=x2./z;
pp3=x3./z;
pp=[pp1;pp2;pp3];
end

function [g,L,LL,Gm2,pp]=SDR_third_alternative(A,B,C) 
n=size(A,1);m=3;
I=eye(n);O=0*I;Ai=inv(A);Ai2=Ai*Ai;Ai3=Ai2*Ai;BB=B*B'; 
T6=[Ai*BB*Ai';O;Ai*BB*Ai3'+Ai3*BB*Ai'-Ai2*BB*Ai2';O]; 
T7=[O,Ai2*BB*Ai2',O,O;O,O,Ai3*BB*Ai3',O;O,O,O,O;O,O,O,O]; 
T4=[Ai;-Ai2;Ai3;O];

cvx_begin sdp; 
cvx_precision best
variable g(1); 
variable S1(n,n); 
variable S2(n,n); 
variable S3(n,n); 
variable S4(n,n); 
variable G12(n,n); 
variable G13(n,n); 
variable G23(n,n);
variable P(n,n);

minimize(g) 
subject to 
S=[S1;S2;S3+P;S4]; 
G=[O,G12-P,G13,P;G12'-P',O,G23,O;G13',G23',O,O;P',O,O,O]; 
S1+S1'>=O;
S2+S2'>=O; 
S3+S3'>=O;
S4+S4'>=O;
G12+G12'>=0;
G13+G13'>=0;
G23+G23'>=0;

L=[g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G];
L>=0; 
cvx_end

S1=full(S1);S2=full(S2);S3=full(S3);S4=full(S4);
G12=full(G12);G13=full(G13);G23=full(G23);
P=full(P);
S=[S1;S2;S3+P;S4]; 
G=[O,G12-P,G13,P;G12'-P',O,G23,O;G13',G23',O,O;P',O,O,O]; 
L=full(L);
LL=L;
L=L(:,1:1+m*n);
L=L(1:1+m*n,:);
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);x3=x(2*n+1:3*n);
z=([Ai;-Ai2;Ai3]'*x-C'); 
pp1=x1./z;
pp2=x2./z;
pp3=x3./z;
pp=[pp1;pp2;pp3];
[t,tt]=eig(LL);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((Gm2-g)/Gm2);
end

function [g,L,LL,Gm2,pp]=SDR_third_alternative_sedumi(A,B,C) 
n=size(A,1);m=3;
I=eye(n);O=0*I;Ai=inv(A);Ai2=Ai*Ai;Ai3=Ai2*Ai;BB=B*B'; 
T6=[Ai*BB*Ai';O;Ai*BB*Ai3'+Ai3*BB*Ai'-Ai2*BB*Ai2';O]; 
T7=[O,Ai2*BB*Ai2',O,O;O,O,Ai3*BB*Ai3',O;O,O,O,O;O,O,O,O]; 
T4=[Ai;-Ai2;Ai3;O];

cvx_begin sdp; 
cvx_precision best
cvx_solver sedumi
variable g(1); 
variable S1(n,n); 
variable S2(n,n); 
variable S3(n,n); 
variable S4(n,n); 
variable G12(n,n); 
variable G13(n,n); 
variable G23(n,n);
variable P(n,n);

minimize(g) 
subject to 
S=[S1;S2;S3+P;S4]; 
G=[O,G12-P,G13,P;G12'-P',O,G23,O;G13',G23',O,O;P',O,O,O]; 
S1+S1'>=O;
S2+S2'>=O; 
S3+S3'>=O;
S4+S4'>=O;
G12+G12'>=0;
G13+G13'>=0;
G23+G23'>=0;

L=[g,C*(T6+S)';(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'-T7-T7'-G];
L>=0; 
cvx_end

S1=full(S1);S2=full(S2);S3=full(S3);S4=full(S4);
G12=full(G12);G13=full(G13);G23=full(G23);
P=full(P);
S=[S1;S2;S3+P;S4]; 
G=[O,G12-P,G13,P;G12'-P',O,G23,O;G13',G23',O,O;P',O,O,O]; 
L=full(L);
LL=L;
L=L(:,1:1+m*n);
L=L(1:1+m*n,:);
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
x=x(2:end);
x1=x(1:n);x2=x(n+1:2*n);x3=x(2*n+1:3*n);
z=([Ai;-Ai2;Ai3]'*x-C'); 
pp1=x1./z;
pp2=x2./z;
pp3=x3./z;
pp=[pp1;pp2;pp3];
[t,tt]=eig(LL);
[y,i]=min(diag(tt)); 
x=t(:,i);x=x/x(1);
Gm2=-x'*[g,C*(T6)';(T6)*C',-(T6)*T4'-T4*(T6)'-T7-T7']*x+g;
Gm2=abs((Gm2-g)/Gm2);
end

function [ss]=p2s(p)
syms s;
ss=solve(s^3-p(1)*s^2+p(2)*s-p(3)==0,s);
ss=double(ss);
end
