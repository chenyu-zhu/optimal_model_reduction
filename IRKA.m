function [Gm,iter_no,sigma] = IRKA(n,m,G,s)
In=eye(n);
Im=eye(m);
A = G.A; B = G.B; C = G.C; D = G.D;
%Assign initiate value of sigma
sigma = s;
%Number of interations in each random search
N = 0;
%Number of random search
R = 0;

flag=1;

%Comparing sigma and eig(-Am) to find local minimizers
while flag==1 && R<=10
    N=N+1;
    for i=1:m 
        Vm(:,i)=(sigma(i)*In-A)\B; 
        Wm(:,i)=(conj(sigma(i))*In-A')\C';
    end
    
    Am=(Wm'*Vm)\(Wm'*A*Vm);
    neg_eigAm=eig(-Am);
    neg_eigAm=sort(neg_eigAm);
    sigma=sort(sigma);
    
    E=sigma-neg_eigAm;
    %Check if it is the local min.
    if abs(sum(E)/sum(sigma))<0.005
        flag=0;
    else
        sigma=neg_eigAm;
    end
    
    if N>=10000
        R=R+1;
        sigma=abs(randn(m,1));
    end
end

%Construct the local minimizer, Gm
Am = Am; 
Bm=(Wm'*Vm)\(Wm'*B); 
Cm=C*Vm;
Dm=0;
Gm=ss(Am,Bm,Cm,Dm); 
iter_no = N;