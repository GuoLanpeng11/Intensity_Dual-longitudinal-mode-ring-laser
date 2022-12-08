function Zout=Z1(x)
syms v u k n;
lamda0=632.8e-9;
lamda1=(3e8*lamda0)/(3e8+(875e6)*lamda0);
k=2*pi/lamda1;
u=1144e6/k;
Y=exp(-(v/u)^2).*((x+1i*k*v).^(-1));
Zout=(1i*k*pi^(-1/2)).*integral(matlabFunction(Y),-inf,inf,'ArrayValued',true);
end

