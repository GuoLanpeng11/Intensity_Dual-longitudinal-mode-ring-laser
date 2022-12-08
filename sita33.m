function asita33=sita33
omiga0=3e8/632.8e-9;
lamda0=632.8e-9;
omiga1=omiga0+876e6;
lamda1=3e8/omiga1;
nu=omiga0;
ra=15.5e6;
rb=41e6;
r=128e6;
B=1.2;
Q=omiga0/1200e6;
nu1=-2000e6+(omiga0+omiga1)/2:20e6:2000e6+(omiga0+omiga1)/2;
nu2=nu1+15e6;
nu3=nu1+640e6;
nu4=nu3+15e6;
Ra=ra.*ones(1,201);
Rb=rb.*ones(1,201);
R=r.*ones(1,201);
k1=2*pi/lamda0;
u1=1200e6/k1;
k2=2*pi/lamda1;
u2=1144e6/k2;




av11=complex(R,omiga0-nu3);
av12=ra;
av13=complex(R,omiga0-nu3);
av21=av11;
av22=av12;
av23=complex(R,nu3-omiga0);
av31=av11;
av32=rb;
av33=av23;
av41=av11;
av42=av32;
av43=av13;
aT12=0;%-1./av12.*(-2/(k1*u1).*(av11./(k1*u1).*Z(av11)-1i));
aT21=1./av22.*(Z(av21)+Z(av23))./(av21+av23);
aT31=1./av32.*(Z(av31)+Z(av33))./(av31+av33);
aT42=0;%-1./av42.*(-2/(k1*u1).*(av41./(k1*u1).*Z(av41)-1i));
i1sita3333=imag(aT12+aT21+aT31+aT42);



bv11=complex(R,omiga1-nu3);
bv12=ra;
bv13=complex(R,omiga1-nu3);
bv21=bv11;
bv22=bv12;
bv23=complex(R,nu3-omiga1);
bv31=bv11;
bv32=rb;
bv33=bv23;
bv41=bv11;
bv42=bv32;
bv43=bv13;
bT12=0;%-1./bv12.*(-2/(k2*u2).*(bv11./(k2*u2).*Z1(bv11)-1i));
bT21=1./bv22.*(Z1(bv21)+Z1(bv23))./(bv21+bv23);
bT31=1./bv32.*(Z1(bv31)+Z1(bv33))./(bv31+bv33);
bT42=0;%-1./bv42.*(-2/(k2*u2).*(bv41./(k2*u2).*Z1(bv41)-1i));
i2sita3333=imag(bT12+bT21+bT31+bT42);

asita33=0.25*ra*rb*nu*B/(Q*imag(Z(r))).*(0.53.*i1sita3333+0.47.*i2sita3333);
end