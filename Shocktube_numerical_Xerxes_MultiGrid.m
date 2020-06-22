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
pR=10;  
rhoR=8;
tF=0.5;
N=[100 200 300]; %Array defining different grid sizes in spatial points
for G=1:length(N)
    % Some constants
    gamma=1.4;
    Courant=0.95; %Courant Number = 0.95<1 to account for non-linearity effects
    x=linspace(-2,2,N(G));
    dx = x(2)-x(1);
    %% Initial Flow Variables(Pressure,Density, Velocity,Energy,Enthalpy)
    p1 = 1;
    rho1 = 1;
    u1 = 0;
    e1 = p1/((gamma-1)*rho1) + (u1^2)/2;
    h1 = e1 + p1/rho1;
    
    p2 = p1*pR;
    rho2 = rho1*rhoR;
    u2=0;
    e2 = p2/((gamma-1)*rho2) + u2^2/2;
    h2 = e2 + p2/rho2;
    %% Apply Initial Flow Variables across x domain (-2<x<2)
    p = zeros(1,N(G));
    p(1:floor((N(G)+1)/2)) = p2; %set initial values on the left of diaphragm
    p(ceil((N(G)+1)/2):N(G)) = p1; %set initial values on the right of diaphragm
    
    rho = zeros(1,N(G));
    rho(1:floor((N(G)+1)/2)) = rho2;
    rho(ceil((N(G)+1)/2):N(G)) = rho1;
    
    u = zeros(1,N(G));
    u(1:floor((N(G)+1)/2)) = u2;
    u(ceil((N(G)+1)/2):N(G)) = u1;
    
    e = zeros(1,N(G));
    e(1:floor((N(G)+1)/2)) = e2;
    e(ceil((N(G)+1)/2):N(G)) = e1;
    
    h = zeros(1,N(G));
    h(1:floor((N(G)+1)/2)) = h2;
    h(ceil((N(G)+1)/2):N(G)) = h1;
    %% Flux-Vector Splitting Scheme
    %Populate U and F matrices
    U=zeros(3,N(G));
    F=zeros(3,N(G));
    U(1,:) = rho;
    U(2,:) = rho.*u;
    U(3,:) = rho.*e;
    F(1,:) = rho.*u;
    F(2,:) = rho.*u.^2+p;
    F(3,:) = rho.*u.*h;
    c = sqrt(gamma.*p./rho);
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
    
    tcount=0; %Counter to track number of iterations performed
    t=0;
    tic 
 while t<tF
     %Calculate maximum time step that satisfies stability condition
     a = abs(u);
     lambda = max(a+c);
     dt = Courant*dx/lambda;
     
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
     U=Unew;    %Update U array with Unew values
     %Reevaluate all the parameters for plotting
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
     
     t=t+dt;  %time march forward by the dt found at this iteration
     tcount=tcount+1;
 end
 disp(tcount) %Display the number of iterations
 %% Collect data from each iteration (new arrays created for each iteration 
 %  due to changing grid size affecting array size)
 if G==1
     x100=x;
     rho100= U(1,:);
     u100=U(2,:)./U(1,:);
     e100=U(3,:)./U(1,:);
     h100= e + p./rho;
     p100=(gamma-1)*rho.*(e-0.5*(u.^2));
     M100= u./c;
     s100= log(p./rho.^gamma);
 elseif G==2
     x200=x;
     rho200= U(1,:);
     u200=U(2,:)./U(1,:);
     e200=U(3,:)./U(1,:);
     h200= e + p./rho;
     p200=(gamma-1)*rho.*(e-0.5*(u.^2));
     M200 = u./c;
     s200= log(p./rho.^gamma);
 else
     x300=x;
     rho300 = U(1,:);
     u300=U(2,:)./U(1,:);
     e300=U(3,:)./U(1,:);
     h300= e + p./rho;
     p300=(gamma-1)*rho.*(e-0.5*(u.^2));
     M300= u./c;
     s300= log(p./rho.^gamma);
 end

end
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
%% Plot numerical vs analytical solutions for each parameter for 3 grid sizes
figure(1)   %Velocity Graph
hold on
plot(x100,u100,'color',[0, 0.4470, 0.7410],'LineWidth',1.25);
plot(x200,u200,'color',[0.9290, 0.6940, 0.1250],'LineWidth',1.25);
plot(x300,u300,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.25);
plot(xAn,uAn,'k:','LineWidth',2.0);
set(gca,'FontSize',13)
grid minor
title('Velocity','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Velocity','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off','FontSize',10)
legend('100 Grid Points','200 Grid Points','300 Grid Points','Analytical');

figure(2)   %Pressure Graph
hold on
plot(x100,p100,'color',[0, 0.4470, 0.7410],'LineWidth',1.25);
plot(x200,p200,'color',[0.9290, 0.6940, 0.1250],'LineWidth',1.25);
plot(x300,p300,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.25);
plot(xAn,pAn,'k:','LineWidth',2.0);
set(gca,'FontSize',13)
grid minor
title('Pressure','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Pressure','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off','FontSize',10)
legend('100 Grid Points','200 Grid Points','300 Grid Points','Analytical');

figure(3)   %Density Graph
hold on
plot(x100,rho100,'color',[0, 0.4470, 0.7410],'LineWidth',1.25);
plot(x200,rho200,'color',[0.9290, 0.6940, 0.1250],'LineWidth',1.25);
plot(x300,rho300,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.25);
plot(xAn,rhoAn,'k:','LineWidth',2.0);
set(gca,'FontSize',13)
grid minor
title('Density','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Density','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off','FontSize',10)
legend('100 Grid Points','200 Grid Points','300 Grid Points','Analytical');

figure(4)   %Local Mach Number Graph
hold on
plot(x100,M100,'color',[0, 0.4470, 0.7410],'LineWidth',1.25);
plot(x200,M200,'color',[0.9290, 0.6940, 0.1250],'LineWidth',1.25);
plot(x300,M300,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.25);
plot(xAn,MAn,'k:','LineWidth',2.0);
set(gca,'FontSize',13)
grid minor
title('Local Mach Number','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Local Mach Number','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off','FontSize',10)
legend('100 Grid Points','200 Grid Points','300 Grid Points','Analytical');

figure(5)   %Entropy Graph
hold on
plot(x100,s100,'color',[0, 0.4470, 0.7410],'LineWidth',1.25);
plot(x200,s200,'color',[0.9290, 0.6940, 0.1250],'LineWidth',1.25);
plot(x300,s300,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.25);
plot(xAn,sAn,'k:','LineWidth',2.0);
set(gca,'FontSize',13)
grid minor
title('Entropy/Cv','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('Entropy/Cv','interpreter','latex');
set(legend,'Interpreter','latex','location','west','box','off','FontSize',10)
legend('100 Grid Points','200 Grid Points','300 Grid Points','Analytical');