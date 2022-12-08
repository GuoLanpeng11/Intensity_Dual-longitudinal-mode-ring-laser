function Zout=Z(x)
syms v;
lamda=632.8e-9;
k=2*pi/lamda;
u=1200e6/k;
Y=exp(-(v/u)^2).*((x+1i*k*v).^(-1));
Zout=(1i*k*pi^(-1/2)).*integral(matlabFunction(Y),-inf,inf,'ArrayValued',true);
end

