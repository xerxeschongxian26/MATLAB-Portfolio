clear
clc
%This script was written to implement the semi-Lagrangian numerical method
%to solve the Advection Equation. 
%This script is specifically for Q10 of the Coursework 1 Handout.
%Written by Mr Xerxes Chong Xian
%CID:01389744
%--------------------------------------------------------------------------
Time =1;
%Constants--------------------
n=101; %Number of grid points
DeltaX=(1-0)/(n-1); %Range of x from 0 to 1 divided by number of intervals
a=0.5;
lambda=[0.5 1.0 1.5 5.0]; % Array of Courant Numbers
x=0:1/(n-1):1; %Discretise x over 101 grid points for 0<x<1
%% Begin semi-Lagrangian Numerical Method
for Q=1:length(lambda)
    DeltaT=lambda(Q)*DeltaX/a;
    tarray=0:DeltaT:Time;                %Discretise an array of time value given delta t
    u=zeros([length(tarray),length(x)]); %Create empty matrix of u values
    for r=1:length(x)                    %Populate u matrix with initial conditions t=0
        if x(r)<=0.5
            u(1,r)=0;
        else
            u(1,r)=1;
        end
    end
    u(1,length(x))=u(1,1);          %Enforce periodic condition
    p=2;                            %Start while loop at the 2nd grid point in time value to enter the while loop
    while p<=length(tarray)
        for i=1:length(x)
            k=floor(i-lambda(Q)+1); %Calculate the value of k based on the inequality relationship
            if k>1
                u(p,i)=-(i-k)*u(p-1,k-1) + (i-k+1)*u(p-1,k)-lambda(Q)*(u(p-1,k)-u(p-1,k-1));
            else
                k=k+100;            %This accounts for periodic conditions when k<=1 due to presence of u_(k-1) terms 
                u(p,i)=-(i-(k-100))*u(p-1,k-1) + (i-(k-100)+1)*u(p-1,k)-lambda(Q)*(u(p-1,k)-u(p-1,k-1));
            end
        end
        u(p,length(x))=u(p,1);      %Enforcing Periodic conditions
        p=p+1;
    end
U(Q,:)=u(length(tarray),:);
end
%% Plot numerical solutions with different Courant Numbers as per Q10
hold on
ylim([0 1.1]);
grid minor
plot(x,U(1,:),'r','color',[0, 0.4470, 0.7410],'LineWidth',1.5)
plot(x,U(2,:),'b--')
plot(x,U(3,:),'g','color',[0.9290, 0.6940, 0.1250],'LineWidth',1.5)
plot(x,U(4,:),'k--','LineWidth',0.75)
legend('0.5','1.0','1.5','5.0');
set(gca,'FontSize',13)
set(legend,'Interpreter','latex','location','northeast','box','off')
title(legend,'Courant Number, $\lambda$')
title('$u(x,t)$ vs $x$','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$u(x,t)$','interpreter','latex');