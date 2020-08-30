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
% numerator=[41 50 140];
% denominator=[1 11 111 110 100];

%system 4
% numerator=[-2.9239 -39.5525 -97.5270 -147.1508];
% denominator=[1 11.9584 43.9119 73.6759 44.3821];

%system 5
% numerator=[-1.2805 -6.2266 -12.8095 -9.3373];
% denominator=[1 3.1855 8.9263 12.2936 3.1987];

%system 6
numerator=[-1.3369 -4.8341 -47.5819 -42.7285];
denominator=[1 17.0728 84.9908 122.4400 59.9309];


sys = tf(numerator,denominator); 
[A,B,C,D] = tf2ss(numerator,denominator); 
G=ss(A,B,C,D);

n=size(A);m=2;
In=eye(n);Im=eye(m);


%Construct new A,B,C,D
O=zeros(n);
I=eye(n);
A11=[A 2*A;O A];
A12=[I O;I O];
A21=[O O;I O];
A22=[A O;O A];

O=zeros(size(B));
B1=[B O;B O];
B2=[O B;O O];

O=zeros(size(C));
C1=[C O;O O];
C2=[O O;O C];

A_new=[A11,A12;A21,A22];
B_new=[B1;B2];
C_new=[C1,C2];
D_new=[zeros(size(C,1),size(B,2)) zeros(size(C,1),size(B,2));
   zeros(size(C,1),size(B,2)) zeros(size(C,1),size(B,2))];

I2N=eye(2*n);
O2N=zeros(2*n);
F1=[I2N O2N;O2N O2N];
F2=[O2N O2N;O2N I2N];

syms s1 s2;

%calculate solutions of p(s1,s2)=0 and p(s2,s1)=0 
pequation(1)=det(D_new-C_new*inv(A_new-s1*F1-s2*F2)*B_new);%Schur complement
pequation(2)=det(D_new-C_new*inv(A_new-s2*F1-s1*F2)*B_new);

[ss1,ss2]=solve(pequation,[s1,s2]);
ss1=double(ss1);
ss2=double(ss2);

for k=1:size(ss1,1)
    if abs(imag(ss1(k)))<0.00001
        ss1(k)=real(ss1(k));
    end
    if abs(imag(ss2(k)))<0.00001
        ss2(k)=real(ss2(k));
    end
end

%remove repeated pairs
for i=1:size(ss1,1)
    for k=i+1:size(ss2,1)
        if ss1(i)==ss2(k) && ss1(k)==ss2(i)
            ss1(k)=0;
            ss2(k)=0;
        end
    end
end

fixedpoints1=[];
fixedpoints2=[];

%remove s1=s2
for k=1:size(ss1,1)
    if ss1(k)~=ss2(k)
        fixedpoints1=[fixedpoints1;ss1(k)];
        fixedpoints2=[fixedpoints2;ss2(k)];
    end
end

%remove negative pairs
positive_reals=true(size(fixedpoints1,1),1);
for i=1:size(fixedpoints1,1)
    if (real(fixedpoints1(i))<0||real(fixedpoints2(i))<0)
        positive_reals(i)=false;
    end
end
fixedpoints1=fixedpoints1(positive_reals);
fixedpoints2=fixedpoints2(positive_reals);

%s1,s2 should be real or in conjugate pairs
conjpairs=true(size(fixedpoints1,1),1);
for i=1:size(fixedpoints1,1)
    if fixedpoints1(i)~=conj(fixedpoints2(i))
        positive_reals(i)=false;
    end
end
fixedpoints1=fixedpoints1(conjpairs);
fixedpoints2=fixedpoints2(conjpairs);

%remove unstable pairs
%Generate the corresponding Gm
stable=true(size(fixedpoints1,1),1);
for i=1:size(fixedpoints1,1)
    Vm=[(fixedpoints1(i)*In-A)\B,(fixedpoints1(i)*In-A)\((fixedpoints2(i)*In-A)\B)];
    Wm=[(fixedpoints1(i)*In-A)'\C',(fixedpoints1(i)*In-A)'\((fixedpoints2(i)*In-A)'\C')];

    Am=(Wm'*Vm)\(Wm'*A*Vm);
    eig_Am=eig(Am);
    for k=1:size(eig_Am)
        if eig_Am(k)>=0
            stable(i)=false;
        end
    end
end
fixedpoints1=fixedpoints1(stable);
fixedpoints2=fixedpoints2(stable);

for i=1:size(fixedpoints1,1)
    Gm= Gm_generate(n,m,G,[fixedpoints1(i),fixedpoints2(i)]);
    E(i) = Error_calculate(n,m,G,Gm,[fixedpoints1(i),fixedpoints2(i)]);
end


function Gm = Gm_generate(n,m,G,sigma)
In=eye(n);Im=eye(m); 
A = G.A;B = G.B;C = G.C;

Vm=[(sigma(1)*In-A)\B,(sigma(1)*In-A)\((sigma(2)*In-A)\B)];
Wm=[(sigma(1)*In-A)'\C',(sigma(1)*In-A)'\((sigma(2)*In-A)'\C')];

Am=(Wm'*Vm)\(Wm'*A*Vm);
Bm=(Wm'*Vm)\(Wm'*B); 
Cm=C*Vm;
Dm=0;
Gm=ss(Am,Bm,Cm,Dm);
end

function E = Error_calculate(n,m,G,Gm,s) 
In=eye(n);Im=eye(m);
A = G.A; B = G.B; C = G.C; D = G.D;
Am = Gm.A; Bm = Gm.B; Cm = Gm.C; Dm = Gm.D; 
P=lyap(A,B*B');
Pm=lyap(Am,Bm*Bm');
E = (C*P*C'-Cm*Pm*Cm')/(C*P*C');
E=sqrt(E);
end
