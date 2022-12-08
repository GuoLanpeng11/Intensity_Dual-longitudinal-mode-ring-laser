function asita32=sita32
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
av13=complex(R,omiga0-nu2);
av21=av11;
av22=av12;
av23=complex(R,nu2-omiga0);
av31=av11;
av32=rb;
av33=av23;
av41=av11;
av42=av32;
av43=av13;
aT11=1/av12.*(Z(av11)+Z(av13))./(av11+av13);
aT22=0;%-1/av22.*((Z(av21)-Z(av23))./(av21-av23));
aT32=0;%-1/av32.*((Z(av31)-Z(av33))./(av31-av33));
aT41=1/av42.*(Z(av41)+Z(av43))./(av41+av43);
i1sita3322=aT11+aT22+aT32+aT41;


bv11=complex(R,omiga0-nu3);
bv12=complex(Ra,nu2-nu3);
bv13=complex(R,omiga0-nu3);
bv21=bv11;
bv22=bv12;
bv23=complex(R,nu2-omiga0);
bv31=bv11;
bv32=complex(Rb,nu2-nu3);
bv33=bv23;
bv41=bv11;
bv42=bv32;
bv43=bv13;
bT13=0.5.*1./(bv11-bv12/2).*((-2/(k1*u1).*(bv11./(k1*u1).*Z(bv11)-1i))-(Z(bv13)-Z(bv12/2))./(bv13-bv12/2));
bT23=0.5.*1./(bv21-bv22/2).*((Z(bv23)-Z(bv21))./(bv23-bv21)-(Z(bv23)-Z(bv22/2))./(bv23-bv22/2));
bT33=0.5.*1./(bv31-bv32/2).*((Z(bv33)-Z(bv31))./(bv33-bv31)-(Z(bv33)-Z(bv32/2))./(bv33-bv32/2));
bT43=0.5.*1./(bv41-bv42/2).*((-2/(k1*u1).*(bv41./(k1*u1).*Z(bv41)-1i))-(Z(bv43)-Z(bv42/2))./(bv43-bv42/2));
i1sita3223=0;%bT13+bT23+bT33+bT43;



cv11=complex(R,omiga1-nu3);
cv13=complex(R,omiga1-nu2);
cv21=cv11;
cv23=complex(R,nu2-omiga1);
cv31=cv11;
cv33=cv23;
cv41=cv11;
cv43=cv13;
cT11=1/ra.*(Z1(cv11)+Z1(cv13))./(cv11+cv13);
cT22=0;%-1/ra.*((Z1(cv21)-Z1(cv23))./(cv21-cv23));
cT32=0;%-1/rb.*((Z1(cv31)-Z1(cv33))./(cv21-cv23));
cT41=1/rb.*(Z1(cv41)+Z1(cv43))./(cv41+cv43);
i2sita3322=cT11+cT22+cT32+cT41;



dv11=complex(R,omiga1-nu3);
dv12=complex(Ra,nu2-nu3);
dv13=complex(R,omiga1-nu3);
dv21=dv11;
dv22=dv12;
dv23=complex(R,nu2-omiga1);
dv31=dv11;
dv32=complex(Rb,nu2-nu1);
dv33=dv23;
dv41=dv11;
dv42=dv32;
dv43=dv13;
dT13=0.5.*1./(dv11-dv12/2).*((-2/(k2*u2).*(dv11./(k2*u2).*Z1(dv11)-1i))-(Z1(dv13)-Z1(dv12/2))./(dv13-dv12/2));
dT23=0.5.*1./(dv21-dv22/2).*((Z1(dv23)-Z1(dv21))./(dv23-dv21)-(Z1(dv23)-Z1(dv22/2))./(dv23-dv22/2));
dT33=0.5.*1./(dv31-dv32/2).*((Z1(dv33)-Z1(dv31))./(dv33-dv31)-(Z1(dv33)-Z1(dv32/2))./(dv33-dv32/2));
dT43=0.5.*1./(dv41-dv42/2).*((-2/(k2*u2).*(dv41./(k2*u2).*Z1(dv41)-1i))-(Z1(dv43)-Z1(dv42/2))./(dv43-dv42/2));
i2sita3223=0;%dT13+dT23+dT33+dT43;

i1sita32=imag(i1sita3322+i1sita3223);
i2sita32=imag(i2sita3322+i2sita3223);
asita32=0.25*ra*rb*nu*B/(Q*imag(Z(r))).*(0.53.*i1sita32+0.47.*i2sita32);

end