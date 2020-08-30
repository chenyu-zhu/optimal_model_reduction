 %Demonstration of Semi-Definite Relaxation Method for m = 1
clear;

%input G(s)
%system 1
% A=[0 0 0 -150;
% 1 0 0 -245 ;
% 0 1 0 -113 ;
% 0 0 1 -19];
% B=[4;1;0;0];
% C=[0 0 0 1];
% D=0;

%system 2
% numerator=[1 15 50];
% denominator=[1 5 33 79 50];

%system 3
% A=[ -0.005 ,-0.99  ;
%     -0.99 ,-5000 ];
% B=[1;100];
% C=[1,100];
% D=0;

%system 4
% numerator=[-1.986,19.17,-0.1606];
% denominator=[1,4.857,14.08,23.02];

%system 5
% numerator=[-0.3556,9.402,5.149,-6.527];
% denominator=[1,9.486,31.24,32.85,5.857];

%system 6
numerator=[0.6355,-5.769,1.119];
denominator=[1,3.696,5.037,1.557];

sys = tf(numerator,denominator); 
[A,B,C,D] = tf2ss(numerator,denominator); 

G=ss(A,B,C,D);

n=size(A);m=1;
In=eye(n);Im=eye(m);

%Implement Semi-definite programming
% [g,L,S,pp] = SDR_first(A,B,C);
[g,L,S,pp] = SDR_first_modified(A,B,C);
[Gm,iter_no,sigma] = IRKA(n,m,G,pp(1));


function [g,L,S,pp]=SDR_first(A,B,C) 
n=size(A,1);I=eye(n);O=0*I;
Ai=inv(A);
T6=[Ai*B*B'*Ai'];
T4=[Ai];
cvx_begin sdp
cvx_precision best
variable g(1);
variable S(n,n); 
minimize(g)
subject to 
S+S'>=0;
[g,C*(T6+S)';
(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)']>=0; 
cvx_end

S=full(S); 
L=[g,C*(T6+S)';
(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'];
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); %the smallest eigenvalue
ne=t(:,i);ne=ne/ne(1); %corresoonding normalized eigenvector
x=ne(2:end)';C1=x(1:n); 
x=x';
z = T4'*x-C';
pp = x./z;
end

function [g,L,S,pp]=SDR_first_modified(A,B,C) 
n=size(A,1);I=eye(n);O=0*I;
Ai=inv(A);
T6=[Ai*B*B'*Ai'];
T4=[Ai];
cvx_begin sdp
cvx_precision best
variable g(1);
variable S(n,n); 
minimize(g)
subject to 
S+S'==0;
[g,C*(T6+S)';
(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)']>=0; 
cvx_end

S=full(S); 
L=[g,C*(T6+S)';
(T6+S)*C',-(T6+S)*T4'-T4*(T6+S)'];
%Construct vector pp
[t,tt]=eig(L);
[y,i]=min(diag(tt)); %the smallest eigenvalue
ne=t(:,i);ne=ne/ne(1); %corresoonding normalized eigenvector
x=ne(2:end)';C1=x(1:n); 
x=x';
z = T4'*x-C';
pp = x./z;
end
