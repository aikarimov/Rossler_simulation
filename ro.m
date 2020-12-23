function Udot = ro(U, Ydata)
% Udot = ro(U, Ydata)
% a = 0.2;
% b = 0.2;
% c = 5.7;

Udot = zeros(3,1);
x = U(1); y = U(2); z = U(3);
ahat = U(4); bhat = U(5); chat = U(6);

xm = Ydata(1);
ym = Ydata(2);
zm = Ydata(3);

e = U(1:3) - Ydata;
u1 = e(1) + (zm - 1)*e(3);
u2 = (1 + ahat)*e(2);
u3 = (1 + x - chat)*e(3);

Udot(1) = -y - z - u1;
Udot(2) = x + ahat*y - u2;
Udot(3) = bhat + z*(x - chat) - u3;

Udot(4) = - ym*e(2);
Udot(5) = -e(3);
Udot(6) = zm*e(3);
end

