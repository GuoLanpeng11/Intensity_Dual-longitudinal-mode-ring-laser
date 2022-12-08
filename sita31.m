function asita31=sita31
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
av13=complex(R,omiga0-nu1);
av21=av11;
av22=av12;
av23=complex(R,nu1-omiga0);
av31=av11;
av32=rb;
av33=av23;
av41=av11;
av42=av32;
av43=av13;
aT12=0;%-1./av12.*(Z(av11)-Z(av13))./(av11-av13);
aT21=1./av22.*(Z(av21)+Z(av23))./(av21+av23);
aT31=1./av32.*(Z(av31)+Z(av33))./(av31+av33);
aT42=0;%-1./av42.*(Z(av41)-Z(av43))./(av41-av43);
i1sita3311=aT12+aT21+aT31+aT42;



bv11=complex(R,omiga0-nu3);
bv12=complex(Ra,nu1-nu3);
bv13=complex(R,omiga0-nu3);
bv21=bv11;
bv22=bv12;
bv23=complex(R,nu1-omiga0);
bv31=bv11;
bv32=complex(Rb,nu1-nu3);
bv33=bv23;
bv41=bv11;
bv42=bv32;
bv43=bv13;
bT12=0;%-1./bv12.*(-2/(k1*u1).*(bv11./(k1*u1).*Z(bv11)-1i));
bT21=1./bv22.*(Z(bv21)+Z(bv23))./(bv21+bv23);
bT31=1./bv32.*(Z(bv31)+Z(bv33))./(bv31+bv33);
bT42=0;%-1./bv42.*(-2/(k1*u1).*(bv41./(k1*u1).*Z(bv41)-1i));
i1sita3113=bT12+bT21+bT31+bT42;



cv11=complex(R,omiga1-nu3);
cv12=ra;
cv13=complex(R,omiga1-nu1);
cv21=cv11;
cv22=cv12;
cv23=complex(R,nu1-omiga1);
cv31=cv11;
cv32=rb;
cv33=cv23;
cv41=cv11;
cv42=cv32;
cv43=cv13;
cT12=0;%-1./cv12.*(Z1(cv11)-Z1(cv13))./(cv11-cv13);
cT21=1./cv22.*(Z1(cv21)+Z1(cv23))./(cv21+cv23);
cT31=1./cv32.*(Z1(cv31)+Z1(cv33))./(cv31+cv33);
dT42=0;%-1./cv42.*(Z1(cv41)-Z1(cv43))./(cv41-cv43);
i2sita3311=cT12+cT21+cT31+dT42;



dv11=complex(R,omiga1-nu3);
dv12=complex(Ra,nu1-nu3);
dv13=complex(R,omiga1-nu3);
dv21=dv11;
dv22=dv12;
dv23=complex(R,nu1-omiga1);
dv31=dv11;
dv32=complex(Rb,nu1-nu3);
dv33=dv23;
dv41=dv11;
dv42=dv32;
dv43=dv13;
dT12=0;%-1./dv12.*(-2/(k2*u2).*(dv11./(k2*u2).*Z1(dv11)-1i));
dT21=1./dv22.*(Z1(dv21)+Z1(dv23))./(dv21+dv23);
dT31=1./dv32.*(Z1(dv31)+Z1(dv33))./(dv31+dv33);
dT42=0;%-1./dv42.*(-2/(k2*u2).*(dv41./(k2*u2).*Z1(dv41)-1i));
i2sita3113=dT12+dT21+dT31+dT42;


i1sita31=imag(i1sita3311+i1sita3113);
i2sita31=imag(i2sita3311+i2sita3113);
asita31=0.25*ra*rb*nu*B/(Q*imag(Z(r))).*(0.53.*i1sita31+0.47.*i2sita31);

end