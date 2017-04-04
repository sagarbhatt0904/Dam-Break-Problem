%McCormack's Method
clear all;
clc;
n=100;
for i=1:n
if i<=50
H_old(i)=10;
H(i)=10;
H1(i)=0;
else
H_old(i)=1;
H(i)=1;
H1(i)=0;
end
v_old(i)=0;
v(i)=0;
v1(i)=0;
Q(i)=0;
Q1(i)=0;
Q_old(i)=0;
E_old(i)=0;
E1(i)=0;
end
L=100;dx=1;
dt=0.1;
k=(dt/dx);
g=9.81;
for t=0:dt:4
for i=2:n-1
Q(i)=(H(i)*v(i));
Q_old(i)=(H_old(i)*v_old(i));
E_old(i)=(H_old(i)*(v_old(i)^2))+(0.5*g*(H_old(i)^2));
E1(i)=(H1(i)*(v1(i)^2))+(0.5*g*(H1(i)^2));
H1(i)=H_old(i)-k*(Q_old(i+1)-Q_old(i));
Q1(i)=Q_old(i)-k*(E_old(i+1)-E_old(i));
H(i)=0.5*(H_old(i)+H1(i)-k*(Q1(i)-Q1(i-1)));
Q(i)=0.5*(Q_old(i)+Q1(i)-k*(E1(i)-E1(i-1)));
end
H(1)=H(2);
H(n)=H(n-1);
Q(1)=Q(2);
Q(n)=Q(n-1);
H1(1)=H1(2);
H1(n)=H1(n-1);
Q1(1)=Q1(2);
Q1(n)=Q1(n-1);
v=Q./H;
v_old=Q_old./H_old;
H_old=H;
Q_old=Q;
end
plot(H);
grid on;
xlabel('-------X (m)-------->');
ylabel('-------H (m)-------->');
title('Height vs x by McCormack Method');
