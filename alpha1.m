function aalpha1=alpha1
omiga0=3e8/632.8e-9;
omiga1=omiga0+876e6;

nu1=-2000e6+(omiga0+omiga1)/2:20e6:2000e6+(omiga0+omiga1)/2;
nu2=nu1+15e6;
nu3=nu1+640e6;
nu4=nu3+15e6;
%nu1=-2000e6+(omiga0+omiga1)/2:20e6:2000e6+(omiga0+omiga1)/2;
%nu2=nu1+1e6;
%nu3=nu1+640e6;
%nu4=nu3+1e6;

B=1.2;
Q1=omiga0/1200e6;
Q2=Q1;%Q2=omiga0/1144e6;
ra=15.5e6;
rb=41e6;
r=128e6;
R=r.*ones(1,201); 

x1=complex(R,omiga0-nu1);  
i1alpha1=0.5.*imag(Z(x1))*omiga0*B/(2*Q1*imag(Z(r)));

x2=complex(R,omiga1-nu1); 
i2alpha1=0.5.*imag(Z1(x2))*omiga0*B/(2*Q2*imag(Z1(r)));

aalpha1=0.53.*i1alpha1+0.47.*i2alpha1-omiga0/2*Q1;
%plot(nu1-(omiga0+omiga1)/2,i1alpha3,'linewidth',1.5);
%hold on
%plot(nu1-(omiga0+omiga1)/2,i2alpha3,'linewidth',1.5);
%hold on
%plot(nu1-(omiga0+omiga1)/2,aalpha3,'linewidth',1.5);
end