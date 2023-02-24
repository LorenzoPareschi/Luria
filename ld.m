function p=ld(beta1,beta2,mu,t,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation from 
% Zheng, Q., Math. Biosciences, 1999
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho=beta1/beta2;

m=(mu/beta1)*(exp(beta1*t)-1);
theta=(mu/beta1)*exp(beta1*t);

for j=1:n+1;
    lambda(j)=theta*(j^(-rho)-(j+1)^(-rho));
end

q(1)=-m;
q(2:n+1)=lambda(1:n);


p=zeros(n+1,1);
p(1)=exp(q(1));

for i=2:n;
    pp=0.0;
    for j=0:i-1
       pp=pp+(i-1-j)*q(i-j)*p(j+1);
    end
    p(i)=(1/(i-1))*pp;
end

