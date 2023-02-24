% From "The kinetic theory of mutation rates"
% by L.Pareschi and G.Toscani, Axioms 2023
%
% Simulation of the Barlett formulation

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

% Parameters

epsilon=0.01; % scaling factor
t=6.7;        % final time
beta1=3;      % grow rate of normal cells
beta2=2.8;    % grow rate of mutants cells
mu=10^(-7);   % per-cell per-unit-time mutation rate
n=300;        % number of mutant cells
dt=0.1;       % time step
N=100000;     % number of particles

% Reference solution

lcexact=lc2(beta1,beta2,mu,t,n);

% Scaling of the variables

T=t/epsilon;
beta2=beta2*epsilon;
mu=mu*epsilon;
beta1=beta1*epsilon;
alpha=beta1+mu;

% Initializations

mc=zeros(1,N);
nc=ones(1,N);
md=zeros(1,N);


Ne(1)=1;
Me(1)=0;
Ma(1)=0;
Mc(1)=0;

k=0;

% Barlett evolution
   
for t=dt:dt:T
    
    % Nm <= N is the number of interacting particles
    
    Nm=rounds(dt*N);
    
    mp=randperm(N);
    mi=mp(1:Nm);
    
    eta=poissrnd(beta1*nc(mi));
    
    zeta=poissrnd(mu*nc(mi));
    
    mp=randperm(N);
    mi=mp(1:Nm);
    
    theta=poissrnd(beta2*mc(mi));
    
    mc(mi)=mc(mi)+theta+zeta;   
    nc(mi)=nc(mi)+eta;
    
    if round(k/50)*50==k
        disp(sprintf('time=%3.2f',t*epsilon));
    end
    
    k=k+1;
    
end

disp(sprintf('time=%3.2f',t*epsilon));
[f,mm]=hist(mc(mc<=n),100);
dm=mm(2)-mm(1);
f=f/(dm*N);
plot(mm,f,':o');
hold on;
plot(lcexact,'-k');
xlabel('mutations $v$');
ylabel('$f(v,t)$');

title(sprintf('Barlett formulation: $\\varepsilon=%g$',epsilon));
legend('Kinetic model','Reference solution');
axis([0 n 0 0.014]);
hold off;
drawnow;

figure;
histogram2(mc(mc<=200),nc(mc<=200),100,'FaceColor','flat','ShowEmptyBins','on');
view(0,90);
title(sprintf('Barlett formulation: $\\varepsilon=%g$',epsilon));
xlabel('mutated cells $v$');
ylabel('normal cells $w$');

