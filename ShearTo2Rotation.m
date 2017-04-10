function [Ra,Rt,S2] = ShearTo2Rotation(h, sx, sy)

H = eye(2);
H(1,2) = h;
S1 = diag([sx,sy]);
A1 = H*S1;

T = sy*sy*(h*h+1) + sx*sx;
T1 = T*T - 4*sx*sx*sy*sy;
tx = sqrt((T-sqrt(T1))/2)
ty = sqrt((T+sqrt(T1))/2)

t1 = 2*h*sy*sx;
t2 = T - 2*sx*sx;

th = atan(t1/t2)/2
al = -atan((ty/tx)*tan(th))

Ra = [cos(al),-sin(al);sin(al),cos(al)];
Rt = [cos(th), -sin(th); sin(th), cos(th)];
S2 = diag([tx,ty]);
A2 = Ra*S2*Rt;
end