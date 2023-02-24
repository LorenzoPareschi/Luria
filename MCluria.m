% From "The kinetic theory of mutation rates"
% by L.Pareschi and G.Toscani, Axioms 2023
%
% Simulation of the Luria-Delbrick and Lea-Coulson formulations

close all;
clear all;

% Latex options

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2);

% Model selection

LD=1;         % 1=Lura-Delbruck, 2=Lea-Coulson
if LD==1
    model='Luria-Delbr\"uck case';
else
    model='Lea-Coulson case';
end

% Parameters

epsilon=0.01; % scaling factor
t=6.7;        % final time
beta1=3;      % grow rate of normal cells
mu=10^(-7);   % per-cell per-unit-time mutation rate
n=300;        % number of mutant cells
dt=0.1;       % time step
N=100000;      % number of particles

% Exact solution

n=300;

if LD==1
    beta2=2.5; % grow rate of mutants cells
    ldexact=ld(beta1,beta2,mu,t,n);
else
    beta2=2.8; % grow rate of mutants cells
    ldexact=lc2(beta1,beta2,mu,t,n);
end

% Scaling of the variables

T=t/epsilon;
beta2=beta2*epsilon;
mu=mu*epsilon;
beta1=beta1*epsilon;
alpha=beta1+mu;

% Initialize variables

mc=zeros(1,N);

Nc(1)=1;
Me(1)=0;
Ma(1)=0;
Vc(1)=0;
Mc(1)=0;

k=1;

for t=dt:dt:T
    k=k+1;
    
    % Evolution of normal cells
    
    Nc(k)=Nc(1)*exp(beta1*t);
    
    % Nm <= N is the number of interacting particles
    
    Nm=rounds(dt*N);
    
    mp=randperm(N);
    mi=mp(1:Nm);
    
    eta=poissrnd(mu*Nc(k-1),1,Nm);
    
    if LD==1
        theta=beta2*mc(mi);
    else
        theta=poissrnd(beta2*mc(mi));
    end
   
    mc(mi)=mc(mi)+theta+eta;

    if round(k/50)*50==k
        disp(sprintf('time=%3.2f',t*epsilon));
    end
    
end

[f,mm]=hist(mc(mc<=n),100);
dm=mm(2)-mm(1);
f=f/(dm*N);
plot(mm,f,':o');
hold on;
plot(ldexact,'-k');
xlabel('mutations $v$');
ylabel('$f(v,t)$');
title(sprintf('%s $:\\varepsilon=%g$',model,epsilon));
legend('Kinetic model','Reference solution');
hold off;
if LD==1
    axis([0 n 0 0.018]);
else
    axis([0 n 0 0.014]);
end
drawnow;
