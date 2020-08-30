clear;

%input G(s)

%system 1
% A=[0 0 0 -150;
% 1 0 0 -245 ;
% 0 1 0 -113 ;
% 0 0 1 -19];
% B=[4;1;0;0];
% C=[0 0 0 1];
% 
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

%Construct tilde_H(s)
tilde_H = ss([A,A;0*A,A],[B;2*B],[C,0*C],0); 
zz = tzero(tilde_H);
%Compute the real and positive zeros of tilde_H(s)
i = find((imag(zz) == 0)&(real(zz)) > 0); 
fixed_points = zz(i);

for i=1:size(fixed_points,1)
    %generate the corresponding reduced system
    Gm= Gm_generate(n,m,G,fixed_points(i));
    %Compute the realtive error
    E(i) = Error_calculate(n,m,G,Gm);
end

function Gm = Gm_generate(n,m,G,sigma)
In=eye(n);Im=eye(m); A = G.A;
B = G.B;
C = G.C;
for i=1:m
    Vm(:,i)=inv(sigma(i)*In-A)*B; 
    Wm(:,i)=inv(conj(sigma(i))*In-A')*C';
end
Am=(Wm'*Vm)\(Wm'*A*Vm);
Bm=(Wm'*Vm)\(Wm'*B); 
Cm=C*Vm;
Dm=0;
Gm=ss(Am,Bm,Cm,Dm);
end

function E = Error_calculate(n,m,G,Gm) 
In=eye(n);Im=eye(m);
A = G.A; B = G.B; C = G.C; D = G.D;
Am = Gm.A; Bm = Gm.B; Cm = Gm.C; Dm = Gm.D; 
P=lyap(A,B*B');
Pm=lyap(Am,Bm*Bm');
E = (C*P*C'-Cm*Pm*Cm')/(C*P*C'); %relative error
E=sqrt(E);
end