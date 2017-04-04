%Analytical Solution
L = 100;
dx = 1;
g=9.81;
t = 4;
n=100;
H1=10;
H2 = 1;
Hg = 3.962;
x1=-sqrt(g*H1)*t;
x2= -sqrt(g*Hg)*t;x3=g*t;
for i=1:n
if (i-L/2)<=x1
H(i) = H1;
u(i) = 0;
end
if ((i-L/2)<=x2 & (i-L/2)>x1)
p=(2*sqrt(g*H1)/3)-(i-L/2)/(3*t);
H(i) = p^2/g;
end
if ((i-L/2)<=x3 & (i-L/2)>x2)
H(i) = Hg;
end
if ((i-L/2)> x3)
H(i) = H2;
end
end
plot(H);
ylim([0 11]);
grid on;
xlabel('-------X (m)-------->');
ylabel('-------H (m)-------->');
title('Height vs x Analytical Solution');
