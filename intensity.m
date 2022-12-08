syms x Y;
omiga=3e8/632.8e-9/1e6;%MHz
v1=-2400+omiga:10:2400+omiga;
v2=v1+0.03;
v3=v1+60;
v4=v3+0.03;
ra=15;%MHz
rb=41;%MHz
r=187;%MHz
sita1111=F(v1,v1).*D(r,-0.5.*v1-0.5.*v1+v1);
sita11=0.5*ra*rb.*imag(sita1111+sita1111);
sita1122=F(v2,v2).*D(r,omiga-0.5.*v1+0.5.*v2-v2);
sita1221=0;
sita12=0.5*ra*rb.*imag(sita1122+sita1221);
sita1133=F(v3,v3).*D(r,-0.5.*v1-0.5.*v3+v3);
sita1331=F(v3,v1).*D(r,-0.5.*v3-0.5.*v1+v3);
sita13=0.5*ra*rb.*imag(sita1133+sita1331);
sita1144=F(v4,v4).*D(r,omiga-0.5.*v1+0.5.*v4-v4);
sita1441=0;
sita14=0.5*ra*rb.*imag(sita1144+sita1441);

sita2211=F(v1,v1).*D(r,omiga-0.5.*v2+0.5.*v1-v1);
sita2112=0;
sita21=0.5*ra*rb.*imag(sita2211+sita2112);
sita2222=F(v2,v2).*D(r,-0.5.*v2-0.5.*v2+v2);
sita22=0.5*ra*rb.*imag(sita2222+sita2222);
sita2233=F(v3,v3).*D(r,omiga-0.5.*v2+0.5.*v3-v3);
sita2332=0;
sita23=0.5*ra*rb.*imag(sita2233+sita2332);
sita2244=F(v4,v4).*D(r,-0.5.*v2-0.5.*v4+v4);
sita2442=F(v4,v2).*D(r,-0.5.*v4-0.5.*v2+v4);
sita24=0.5*ra*rb.*imag(sita2244+sita2442);

sita3311=F(v1,v1).*D(r,-0.5.*v3-0.5.*v1+v1);
sita3113=F(v1,v3).*D(r,-0.5.*v1-0.5.*v3+v1);
sita31=0.5*ra*rb.*imag(sita3311+sita3113);
sita3322=F(v2,v2).*D(r,omiga-0.5.*v3+0.5.*v2-v2);
sita3223=0;
sita32=0.5*ra*rb.*imag(sita3322+sita3223);
sita3333=F(v3,v3).*D(r,-0.5.*v3-0.5.*v3+v3);
sita33=0.5*ra*rb.*imag(sita3333+sita3333);
sita3344=F(v4,v4).*D(r,omiga-0.5.*v3+0.5.*v4-v4);
sita3443=0;
sita34=0.5*ra*rb.*imag(sita3344+sita3443);

sita4411=F(v1,v1).*D(r,omiga-0.5.*v4+0.5.*v1-v1);
sita4114=0;
sita41=0.5*ra*rb.*imag(sita4411+sita4114);
sita4422=F(v2,v2).*D(r,-0.5.*v4-0.5.*v2+v2);
sita4224=F(v2,v4).*D(r,-0.5.*v2-0.5.*v4+v2);
sita42=0.5*ra*rb.*imag(sita4422+sita4224);
sita4433=F(v3,v3).*D(r,omiga-0.5.*v4+0.5.*v3-v3);
sita4334=0;
sita43=0.5*ra*rb.*imag(sita4433+sita4334);
sita4444=F(v4,v4).*D(r,-0.5.*v4-0.5.*v4+v4);
sita44=0.5*ra*rb.*imag(sita4444+sita4444);

ku=1100;
kxi1=(omiga-v1)./ku;
Y=@(x)exp(-x.^2);
alpha1=imagZ(kxi1).*(1./imagZ(r));

kxi2=(omiga-v2)./ku;
Y=@(x)exp(-x.^2);
alpha2=imagZ(kxi2).*(1./imagZ(r));

kxi3=(omiga-v3)./ku;
Y=@(x)exp(-x.^2);
alpha3=imagZ(kxi3).*(1./imagZ(r));

kxi4=(omiga-v4)./ku;
Y=@(x)exp(-x.^2);
alpha4=imagZ(kxi4).*(1./imagZ(r));

 
I1=(alpha1.*(sita22.*sita33.*sita44 - sita44.*sita23.*sita32 - sita33.*sita24.*sita42 - sita22.*sita34.*sita43 + sita23.*sita34.*sita42 + sita24.*sita32.*sita43))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha4.*(sita22.*sita33.*sita14 - sita33.*sita12.*sita24 - sita22.*sita13.*sita34 + sita12.*sita23.*sita34 + sita13.*sita24.*sita32 - sita14.*sita23.*sita32))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha3.*(sita22.*sita44.*sita13 - sita44.*sita12.*sita23 - sita22.*sita14.*sita43 + sita12.*sita24.*sita43 - sita13.*sita24.*sita42 + sita14.*sita23.*sita42))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha2.*(sita33.*sita44.*sita12 - sita44.*sita13.*sita32 - sita33.*sita14.*sita42 - sita12.*sita34.*sita43 + sita13.*sita34.*sita42 + sita14.*sita32.*sita43))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41);
%subplot(2,2,1);
plot(v1-omiga,I1,'r','linewidth',1.5);
hold on
%xlabel('失谐频率(MHz)');
%ylabel('光强I_1');
 
I2=(alpha2.*(sita11.*sita33.*sita44 - sita44.*sita13.*sita31 - sita33.*sita14.*sita41 - sita11.*sita34.*sita43 + sita13.*sita34.*sita41 + sita14.*sita31.*sita43))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha4.*(sita11.*sita33.*sita24 - sita33.*sita14.*sita21 - sita11.*sita23.*sita34 + sita13.*sita21.*sita34 - sita13.*sita24.*sita31 + sita14.*sita23.*sita31))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha3.*(sita11.*sita44.*sita23 - sita44.*sita13.*sita21 - sita11.*sita24.*sita43 + sita13.*sita24.*sita41 + sita14.*sita21.*sita43 - sita14.*sita23.*sita41))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha1.*(sita33.*sita44.*sita21 - sita44.*sita23.*sita31 - sita33.*sita24.*sita41 - sita21.*sita34.*sita43 + sita23.*sita34.*sita41 + sita24.*sita31.*sita43))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41);
%subplot(2,2,2);
plot(v1-omiga,I2,'c','linewidth',1.5);
hold on
%xlabel('失谐频率(MHz)');
%ylabel('光强I_2');
 
I3=(alpha3.*(sita11.*sita22.*sita44 - sita44.*sita12.*sita21 - sita22.*sita14.*sita41 - sita11.*sita24.*sita42 + sita12.*sita24.*sita41 + sita14.*sita21.*sita42))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha4.*(sita11.*sita22.*sita34 - sita22.*sita14.*sita31 - sita11.*sita24.*sita32 - sita12.*sita21.*sita34 + sita12.*sita24.*sita31 + sita14.*sita21.*sita32))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha2.*(sita11.*sita44.*sita32 - sita44.*sita12.*sita31 - sita11.*sita34.*sita42 + sita12.*sita34.*sita41 + sita14.*sita31.*sita42 - sita14.*sita32.*sita41))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha1.*(sita22.*sita44.*sita31 - sita44.*sita21.*sita32 - sita22.*sita34.*sita41 + sita21.*sita34.*sita42 - sita24.*sita31.*sita42 + sita24.*sita32.*sita41))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41);
%subplot(2,2,3)
plot(v1-omiga,I3,'b','linewidth',1.5); 
hold on
%xlabel('失谐频率(MHz)');
%ylabel('光强I_3');
 
I4=(alpha4.*(sita11.*sita22.*sita33 - sita33.*sita12.*sita21 - sita22.*sita13.*sita31 - sita11.*sita23.*sita32 + sita12.*sita23.*sita31 + sita13.*sita21.*sita32))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha3.*(sita11.*sita22.*sita43 - sita22.*sita13.*sita41 - sita11.*sita23.*sita42 - sita12.*sita21.*sita43 + sita12.*sita23.*sita41 + sita13.*sita21.*sita42))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha2.*(sita11.*sita33.*sita42 - sita33.*sita12.*sita41 - sita11.*sita32.*sita43 + sita12.*sita31.*sita43 - sita13.*sita31.*sita42 + sita13.*sita32.*sita41))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41) - (alpha1.*(sita22.*sita33.*sita41 - sita33.*sita21.*sita42 - sita22.*sita31.*sita43 + sita21.*sita32.*sita43 + sita23.*sita31.*sita42 - sita23.*sita32.*sita41))./(sita11.*sita22.*sita33.*sita44 - sita33.*sita44.*sita12.*sita21 - sita22.*sita44.*sita13.*sita31 - sita11.*sita44.*sita23.*sita32 - sita22.*sita33.*sita14.*sita41 - sita11.*sita33.*sita24.*sita42 - sita11.*sita22.*sita34.*sita43 + sita44.*sita12.*sita23.*sita31 + sita44.*sita13.*sita21.*sita32 + sita33.*sita12.*sita24.*sita41 + sita33.*sita14.*sita21.*sita42 + sita22.*sita13.*sita34.*sita41 + sita22.*sita14.*sita31.*sita43 + sita11.*sita23.*sita34.*sita42 + sita11.*sita24.*sita32.*sita43 + sita12.*sita21.*sita34.*sita43 - sita12.*sita23.*sita34.*sita41 - sita12.*sita24.*sita31.*sita43 - sita13.*sita21.*sita34.*sita42 + sita13.*sita24.*sita31.*sita42 - sita13.*sita24.*sita32.*sita41 - sita14.*sita21.*sita32.*sita43 - sita14.*sita23.*sita31.*sita42 + sita14.*sita23.*sita32.*sita41);
%subplot(2,2,4);
plot(v1-omiga,I4,'g','linewidth',1.5);
hold on
%xlabel('失谐频率(MHz)');
%ylabel('光强I_4');

%Is=(I1+I3)/2;
%subplot(2,3,5);
%plot(v1-omiga,Is,'c','linewidth',1.5);
%xlabel('失谐频率(MHz)');
%ylabel('强模平均光强(I_1+I_3)/2');

%Iw=(I2+I4)/2;
%subplot(2,3,6);
%plot(v1-omiga,Iw,'r','linewidth',1.5);
%xlabel('失谐频率(MHz)');
%ylabel('弱模平均光强(I_2+I_4)/2');

