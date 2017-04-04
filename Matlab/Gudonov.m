% Gudunov's Method
clear all;
clc;
n=100;
for i=1:n
if i<=50
H_old(i)=10;
H(i)=10;
else
H_old(i)=1;
H(i)=1;
end
v_old(i)=0;
v(i)=0;
Q(i)=0;
Q_old(i)=0;
E_old(i)=0;
end
L=100;
dx=1;
dt=0.1;
k=(dt/dx);
g=9.81;
for t=0:dt:4
for i=2:n-1
Q(i)=(H(i)*v(i));
Q_old(i)=(H_old(i)*v_old(i));
E_old(i)=(H_old(i)*(v_old(i)^2))+(0.5*g*(H_old(i)^2));
a(i)=max(abs(v_old(i)+sqrt(g*H_old(i))),abs(v_old(i-
1)+sqrt(g*H_old(i-1))));
F1(i)=0.5*(Q_old(i)+Q_old(i-1))-0.5*(H_old(i+1)-H_old(i));
F2(i)=0.5*(E_old(i)+E_old(i-1))-0.5*(Q_old(i+1)-Q_old(i));
end
for i=2:n-1
H(i)=H_old(i)-k*(F1(i)-F1(i-1));
Q(i)=Q_old(i)-k*(F2(i)-F2(i-1));
end
H=smooth(H);
Q=smooth(Q);
H(1)=H(2);
H(n)=H(n-1);
Q(1)=Q(2);
Q(n)=Q(n-1);H_old=H;
Q_old=Q;
v_old=v;
v=Q./H;
end
for i=1:n
if H(i)>10
H(i)=10;
end
end
plot(H);
grid on;
xlabel('-------X (m)-------->');
ylabel('-------H (m)-------->');
title('Height vs x by Gudunov');
