clear
clc
%This script was written to implement the semi-Lagrangian numerical method
%to solve the Advection Equation. 
%This script is specifically for Q11 of the Coursework 1 Handout.
%Written by Mr Xerxes Chong Xian
%CID:01389744
%--------------------------------------------------------------------------
%% Collect inputs from user and define constants
Time =1;
n=501;                %Vary the number of grid points here
DeltaX=1/(n-1);
DeltaT=0.1;
tarray=0:DeltaT:Time; %Discretise an array of time value given DeltaT
x=0:DeltaX:1;         %Discretise x over n grid points for 0<x<1
u=zeros([length(tarray),length(x)]); %Create empty matrix of u values
for r=1:length(x)     %Populate u matrix with initial conditions t=0
    if x(r)<=0.5
        u(1,r)=0;
    else
        u(1,r)=1;
    end
end
u(1,length(x))=u(1,1);           %Enforce periodic condition
p=2;                             %Start while loop at the 2nd grid point in time value to enter the while loop
while p<=length(tarray)
    a=10*sin(tarray(p));
    lambda=a*DeltaT/DeltaX;      %Courant-Number
        for i=1:length(x)
            k=floor(i-lambda+1); %Calculate the value of k at each spacial grid point i
            if k>1
                u(p,i)=-(i-k)*u(p-1,k-1) + (i-k+1)*u(p-1,k)-lambda*(u(p-1,k)-u(p-1,k-1));
            else
                k=k+(n-1);
                u(p,i)=-(i-(k-(n-1)))*u(p-1,k-1) + (i-(k-(n-1))+1)*u(p-1,k)-lambda*(u(p-1,k)-u(p-1,k-1));
            end
        end
    p=p+1;
end
%% Plot numerical solutions with varying grid points as per Q11
hold on
grid minor
plot(x,u(length(tarray),:),'color',[0.25, 0.25, 0.25],'LineWidth',1.0)
legend('11','21','51','101','201','251','501');
set(gca,'FontSize',13)
set(legend,'Interpreter','latex','location','northwest','box','off')
title(legend,'Number of Grid Points')
title('$u(x,t)$ vs $x$','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$u(x,t)$','interpreter','latex');