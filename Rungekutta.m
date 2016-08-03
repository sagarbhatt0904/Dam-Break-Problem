% Runge-Kutta Method
clear all;
clc;
n=100;
for i=1:n
if i<=50
H_old(i)=10;
H(i)=10;
H1(i)=10;
H2(i)=10;
H3(i)=10;
H4(i)=10;
else
H_old(i)=1;
H(i)=1;H1(i)=1;
H2(i)=1;
H3(i)=1;
H4(i)=1;
end
v_old(i)=0;
v(i)=0;
Q(i)=0;
Q_old(i)=0;
E_old(i)=0;
Q1(i)=0;
Q2(i)=0;
Q3(i)=0;
Q4(i)=0;
E1(i)=0;
E2(i)=0;
E3(i)=0;
E4(i)=0;
a_old_H(i)=0;
a_H(i)=0;
sigma_old_H(i)=0;
sigma_H(i)=0;
sigma_old_Q(i)=0;
sigma_Q(i)=0;
theta_old_H(i)=0;
theta_H(i)=0;
theta_old_Q(i)=0;
theta_Q(i)=0;
G_old_H(i)=0;
G_H(i)=0;
G_old_Q=0;
G_Q(i)=0;
psi_old_H(i)=0;
psi_H(i)=0;
psi_old_Q(i)=0;
psi_Q(i)=0;
end
L=100;
dx=1;
dt=0.001;
k=(dt/dx);
g=9.81;
for t=0:dt:4
for i=2:n-2
if H_old(i)-H_old(i-1)<=0.0005
a_old_H(i)=H_old(i);
else
a_old_H(i)=((Q_old(i)-Q_old(i-1))/(H_old(i)-H_old(i-1)));
end
if H_old(i+1)-H_old(i)<=0.0005
a_H(i)=H_old(i);
else
a_H(i)=((Q_old(i+1)-Q_old(i))/(H_old(i+1)-H_old(i)));
endif Q_old(i)-Q_old(i-1)<=0.0005
a_old_Q(i)=Q_old(i);
else
a_old_Q(i)=((E_old(i)-E_old(i-1))/(Q_old(i)-Q_old(i-1)));
end
if Q_old(i+1)-Q_old(i)<=0.0005
a_Q(i)=Q_old(i);
else
a_Q(i)=((E_old(i+1)-E_old(i))/(Q_old(i+1)-Q_old(i)));
end
sigma_old_H(i)=sign(a_old_H(i));
sigma_H(i)=sign(a_H(i));
sigma_old_Q(i)=sign(a_old_Q(i));
sigma_Q(i)=sign(a_Q(i));
if H_old(i)-H_old(i-1)==0
theta_old_H(i)=0;
else
theta_old_H(i)=(H_old(i+sigma_old_H(i))-H_old(i-
1+sigma_old_H(i)))/(H_old(i)-H_old(i-1));
end
if H_old(i+1)-H_old(i)==0
theta_H(i)=0;
else
theta_H(i)=(H_old(i+1+sigma_H(i))-
H_old(i+sigma_H(i)))/(H_old(i+1)-H_old(i));
end
if Q_old(i)-Q_old(i-1)==0
theta_old_Q(i)=0;
else
theta_old_Q(i)=(Q_old(i+sigma_old_Q(i))-Q_old(i-
1+sigma_old_Q(i)))/(Q_old(i)-Q_old(i-1));
end
if Q_old(i+1)-Q_old(i)==0
theta_Q(i)=0;
else
theta_Q(i)=(Q_old(i+1+sigma_Q(i))-
Q_old(i+sigma_Q(i)))/(Q_old(i+1)-Q_old(i));
end
G_old_H(i)=(theta_old_H(i)+abs(theta_old_H(i)))/(1+theta_old_H(i));
G_H(i)=(theta_H(i)+abs(theta_H(i)))/(1+theta_H(i));
G_old_Q(i)=(theta_old_Q(i)+abs(theta_old_Q(i)))/(1+theta_old_Q(i));
G_Q(i)=(theta_Q(i)+abs(theta_Q(i)))/(1+theta_Q(i));
psi_old_H(i)=(0.5*G_old_H(i)*(abs(a_old_H(i))+k*(a_old_H(i)^2))-
abs(a_old_H(i)))*(H_old(i)-H_old(i-1));
psi_H(i)=(0.5*G_H(i)*(abs(a_H(i))+k*(a_H(i)^2))-
abs(a_H(i)))*(H_old(i+1)-H_old(i));
psi_old_Q(i)=(0.5*G_old_Q(i)*(abs(a_old_Q(i))+k*(a_old_Q(i)^2))-
abs(a_old_Q(i)))*(Q_old(i)-Q_old(i-1));
psi_Q(i)=(0.5*G_Q(i)*(abs(a_Q(i))+k*(a_Q(i)^2))-
abs(a_Q(i)))*(Q_old(i+1)-Q_old(i));
end
for i=1:nQ(i)=(H(i)*v(i));
Q_old(i)=(H_old(i)*v_old(i));
E_old(i)=(H_old(i)*(v_old(i)^2))+(0.5*g*(H_old(i)^2));
Q1(i)=(H1(i)*v_old(i));
Q2(i)=(H2(i)*v_old(i));
Q3(i)=(H3(i)*v_old(i));
Q4(i)=(H4(i)*v_old(i));
E1(i)=(H1(i)*(v_old(i)^2))+(0.5*g*(H1(i)^2));
E2(i)=(H2(i)*(v_old(i)^2))+(0.5*g*(H2(i)^2));
E3(i)=(H3(i)*(v_old(i)^2))+(0.5*g*(H3(i)^2));
E4(i)=(H4(i)*(v_old(i)^2))+(0.5*g*(H4(i)^2));
end
for i=2:n-1
H1(i)=H_old(i);
H2(i)=H_old(i)-0.25*0.5*k*(Q1(i+1)-Q1(i-1));
H3(i)=H_old(i)-0.33*0.5*k*(Q2(i+1)-Q2(i-1));
H4(i)=H_old(i)-0.5*0.5*k*(Q3(i+1)-Q3(i-1));
H(i)=H_old(i)-0.5*k*(Q4(i+1)-Q4(i-1));
Q1(i)=Q_old(i);
Q2(i)=Q_old(i)-0.25*0.5*k*(E1(i+1)-E1(i-1));
Q3(i)=Q_old(i)-0.33*0.5*k*(E2(i+1)-E2(i-1));
Q4(i)=Q_old(i)-0.5*0.5*k*(E3(i+1)-E3(i-1));
Q(i)=Q_old(i)-0.5*k*(E4(i+1)-E4(i-1));
end
for i=1:n
H(i)=H(i)-0.5*k*(psi_H(i)-psi_old_H(i));
Q(i)=Q(i)-0.5*k*(psi_Q(i)-psi_old_Q(i));
end
H(1)=H(2);
H(n)=H(n-1);
Q(1)=Q(2);
Q(n)=Q(n-1);
v=Q./H;
H_old=H;
Q_old=Q;
v_old=v;
end
H=smooth(H);
Q=smooth(Q);
plot(H);
grid on;
xlabel('-------X (m)-------->');
ylabel('-------H (m)-------->');
title('Height vs x by Runge Kutta Method with TVD');
