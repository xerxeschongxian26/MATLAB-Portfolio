clear
clc
%This script solves analytically the 1D Euler Equations using the
%Flux-Vector Splitting Method by Steger and Warming 
%Written by Mr Xerxes Chong Xian
%CID:01389744
%--------------------------------------------------------------------------
%% Collect inputs from user and define constants
%  pRat = input('Input Pressure Ratio (P_2/P_1): ');
%  rhoRat = input('Input Density Ratio (rho_2/rho_1): ');
%  tFin = input('Input Time Evaluated (t): ');
%  N = input('Input spatial points N: ');
pRat=10;
rhoRat=8;
tFin=0.5;
N=300;
% Define constants & variables
gamma = 1.4;
sig = 0.95;         %Courant Number = 1 for stability
x=linspace(-2,2,N);
dx = x(2)-x(1);

%% Initialise Initial Conditions (Pressure,Density, Velocity, Energy,Enthalpy)
p1 = 1;
rho1 = 1;
u1 = 0;
e1 = p1/((gamma-1)*rho1) + (u1^2)/2;
h1 = e1 + p1/rho1;

p2 = p1*pRat;
rho2 = rho1*rhoRat;
u2=0;
e2 = p2/((gamma-1)*rho2) + u2^2/2;
h2 = e2 + p2/rho2;
%pressure
p = zeros(1,N);
p(1:floor((N+1)/2)) = p2; %set initial values on the left
p(ceil((N+1)/2):N) = p1; %set initial values on the right
%density
rho = zeros(1,N);
rho(1:floor((N+1)/2)) = rho2; %set initial values of rho on the left
rho(ceil((N+1)/2):N) = rho1; %set initial values of rho on the right
%velocity
u = zeros(1,N);
u(1:floor((N+1)/2)) = u2;
u(ceil((N+1)/2):N) = u1;
%energy
e = zeros(1,N);
e(1:floor((N+1)/2)) = e2;
e(ceil((N+1)/2):N) = e1;
%enthalpy
h = zeros(1,N);
h(1:floor((N+1)/2)) = h2;
h(ceil((N+1)/2):N) = h1;
%% Begin Flux-Vector Splitting Scheme
%Generate U and F matrices
U(1,:) = rho;
U(2,:) = rho.*u;
U(3,:) = rho.*e;
F(1,:) = rho.*u;
F(2,:) = rho.*u.^2+p;
F(3,:) = rho.*u.*h;
c = sqrt(gamma.*p./rho); %speed of sound , Varies spatially?
%Obtain the eigenvalues
l1=u-c;
l2=u;
l3=u+c;
%obtain the positive(Lp1,Lp2,Lp3) and negative (Lm1,Lm2,Lm3) eigenvalues
Lp1=0.5*(l1+abs(l1));
Lp2=0.5*(l2+abs(l2));
Lp3=0.5*(l3+abs(l3));
Lm1=0.5*(l1-abs(l1));
Lm2=0.5*(l2-abs(l2));
Lm3=0.5*(l3-abs(l3));

%  a = abs(u);
%  lambda = max(a+c);
%  dt = 0.004; %sig*dx/lambda;
 tcount=1;
 t=0;
 tic
 while t<tFin
     a = abs(u);
     lambda = max(a+c);
     dt = sig*dx/lambda;
     for i=1:length(x)
         Fpi=rho(i)/(2*gamma)*[Lp1(i)+2*(gamma-1)*Lp2(i)+Lp3(i);(u(i)-c(i))*Lp1(i)+2*(gamma-1)*u(i)*Lp2(i)+(u(i)+c(i))*Lp3(i);(h(i)-u(i)*c(i))*Lp1(i)+(gamma-1)*(u(i)^2)*Lp2(i)+(h(i)+u(i)*c(i))*Lp3(i)];
         if i==1 %Boundary Conditions at start
             Fpim1=rho(i)/(2*gamma)*[Lp1(i)+2*(gamma-1)*Lp2(i)+Lp3(i);(u(i)-c(i))*Lp1(i)+2*(gamma-1)*u(i)*Lp2(i)+(u(i)+c(i))*Lp3(i);(h(i)-u(i)*c(i))*Lp1(i)+(gamma-1)*(u(i)^2)*Lp2(i)+(h(i)+u(i)*c(i))*Lp3(i)];
         else
             Fpim1=rho(i-1)/(2*gamma)*[Lp1(i-1)+2*(gamma-1)*Lp2(i-1)+Lp3(i-1);(u(i-1)-c(i-1))*Lp1(i-1)+2*(gamma-1)*u(i-1)*Lp2(i-1)+(u(i-1)+c(i-1))*Lp3(i-1);(h(i-1)-u(i-1)*c(i-1))*Lp1(i-1)+(gamma-1)*(u(i-1)^2)*Lp2(i-1)+(h(i-1)+u(i-1)*c(i-1))*Lp3(i-1)];
         end
         Fmi=rho(i)/(2*gamma)*[Lm1(i)+2*(gamma-1)*Lm2(i)+Lm3(i);(u(i)-c(i))*Lm1(i)+2*(gamma-1)*u(i)*Lm2(i)+(u(i)+c(i))*Lm3(i);(h(i)-u(i)*c(i))*Lm1(i)+(gamma-1)*(u(i)^2)*Lm2(i)+(h(i)+u(i)*c(i))*Lm3(i)];
         if i==length(x) %Boundary Conditions at end
             Fmi1=rho(i)/(2*gamma)*[Lm1(i)+2*(gamma-1)*Lm2(i)+Lm3(i);(u(i)-c(i))*Lm1(i)+2*(gamma-1)*u(i)*Lm2(i)+(u(i)+c(i))*Lm3(i);(h(i)-u(i)*c(i))*Lm1(i)+(gamma-1)*(u(i)^2)*Lm2(i)+(h(i)+u(i)*c(i))*Lm3(i)];
         else
             Fmi1=rho(i+1)/(2*gamma)*[Lm1(i+1)+2*(gamma-1)*Lm2(i+1)+Lm3(i+1);(u(i+1)-c(i+1))*Lm1(i+1)+2*(gamma-1)*u(i+1)*Lm2(i+1)+(u(i+1)+c(i+1))*Lm3(i+1);(h(i+1)-u(i+1)*c(i+1))*Lm1(i+1)+(gamma-1)*(u(i+1)^2)*Lm2(i+1)+(h(i+1)+u(i+1)*c(i+1))*Lm3(i+1)];
         end
         Unew(:,i)=U(:,i)-(dt/dx)*(Fpi-Fpim1 + Fmi1-Fmi);%update the values of U at location i at each time step, U on the left is the same position on the right but at the next time step
     end
     U=Unew; %Update U array with Unew values
     %Reevaluate all the parameters
     rho = U(1,:);
     u=U(2,:)./U(1,:);
     e=U(3,:)./U(1,:);
     h = e + p./rho;
     p=(gamma-1)*rho.*(e-0.5*(u.^2));
     c = sqrt(gamma*p./rho);
     %Obtain the eigenvalues
     l1=u-c;
     l2=u;
     l3=u+c;
     %obtain the positive(Lp1,Lp2,Lp3) and negative (Lm1,Lm2,Lm3) eigenvalues
     Lp1=0.5*(l1+abs(l1));
     Lp2=0.5*(l2+abs(l2));
     Lp3=0.5*(l3+abs(l3));
     Lm1=0.5*(l1-abs(l1));
     Lm2=0.5*(l2-abs(l2));
     Lm3=0.5*(l3-abs(l3));
     
     
     
     tcount=tcount+1;
     t=t+dt; %time march forward by the prescribed dt
 end
 disp(tcount)
%% Generate data for plotting
rho = U(1,:);
u=U(2,:)./U(1,:);
e=U(3,:)./U(1,:);
h = e + p./rho;
p=(gamma-1)*rho.*(e-0.5*(u.^2));
M = u./c;
s = log(p./rho.^gamma);
%% Import analytical data from shock_analytic.dat into an array
shockanalytic= readtable('shock_analytic.dat');
shockanalytic=shockanalytic(:,1:6);  %Appending erroneous last roll
shockanalytic=table2array(shockanalytic);
xAn=shockanalytic(:,1);
rhoAn=shockanalytic(:,2);
pAn=shockanalytic(:,3);
uAn=shockanalytic(:,4);
MAn=shockanalytic(:,5);
sAn=shockanalytic(:,6);
%% Plot numerical vs analytical solutions for each parameter
figure(1)   %Velocity Graph
hold on
plot(x,u,'r','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(xAn,uAn,'r--','LineWidth',2.0);
set(gca,'FontSize',13)
grid minor
title('Velocity','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Velocity','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off')
legend('Numerical','Analytical');


figure(2)   %Pressure Graph
hold on
plot(x,p,'r','color',[0, 0.4470, 0.7410],'LineWidth',1.5)
plot(xAn,pAn,'r--','LineWidth',2.0);
set(gca,'FontSize',13)
grid minor
title('Pressure','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Pressure','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off')
legend('Numerical','Analytical');

figure(3)   %Density Graph
hold on
plot(x,rho,'r','color',[0, 0.4470, 0.7410],'LineWidth',1.5)
plot(xAn,rhoAn,'r--','LineWidth',2.0);
set(gca,'FontSize',13)
grid minor
title('Density','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Density','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off')
legend('Numerical','Analytical');

figure(4)   %Local Mach Number Graph
hold on
plot(x,M,'r','color',[0, 0.4470, 0.7410],'LineWidth',1.5)
plot(xAn,MAn,'r--','LineWidth',2.0);
grid minor
set(gca,'FontSize',13)
title('Local Mach Number','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Local Mach Number','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off')
legend('Numerical','Analytical');

figure(5)   %Entropy Graph
grid minor
hold on
plot(x,s,'r','color',[0, 0.4470, 0.7410],'LineWidth',1.5)
plot(xAn,sAn,'r--','LineWidth',2.0);
set(gca,'FontSize',13)
title('Entropy/CV','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Entropy/CV','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off')
legend('Numerical','Analytical');