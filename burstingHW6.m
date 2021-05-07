
simTime=100;
deltaT=.01;
t=0:deltaT:simTime;

changeTimes=(0);
currentLevels=(50);

I(1:500)=currentLevels; I(501:2000)=0; I(2001:numel(t))=currentLevels;
%Resting potentials
gK=3;gNa=120;
gL=.012;gCa=3.2;gKCa=.02;
EK=-12; ENa=115; EL=10.6;
C=1;
Kd=1;k1=.0275;kc=.02;f=.007;
%state variables
%v=;
%m=;
%n=;
%h=;
%alphas/betas initial
V=0;
alpha_n=.01*((10-V)/(exp((10-V)/10)-1));
beta_n=.125*exp(-V/80);
alpha_m=.1*((25-V)/(exp((25-V)/10)-1));
beta_m=4*exp(-V/18);
alpha_h=.07*exp(-V/20);
beta_h= 1/exp((30-V/10)+1);

n(1)= alpha_n/(alpha_n+beta_n);
m(1)=alpha_m/(alpha_m+beta_m);
h(1)=alpha_h/(alpha_h+beta_h);   

for i=1:numel(t)-1
    %V=0
    alpha_n(i)=.01*((10-V(i))/(exp((10-V(i))/10)-1));
    beta_n(i)=.125*exp(-V(i)/80);
    alpha_m(i)=.1*((25-V(i))/(exp((25-V(i))/10)-1));
    beta_m(i)=4*exp(-V(i)/18);
    alpha_h(i)=.07*exp(-V(i)/20);
    beta_h(i)= 1/(exp((30-V(i))/10)+1);

    %currents
    I_Na=(m(i)^3)*gNa*h(i)*(V(i)-ENa);
    I_Ca=(m(i)^3)*gNa*h(i)*(V(i)-ENa);
    I_K=(n(i)^4)*gK*(V(i)-EK);
    I_L=gL*(V(i)-EL);
       
    %calcium concentration
    c=f*(-k1*I_Ca-kc*c)
    
    %derivatives
    V(i+1)=V(i)+deltaT*I_ion/C;
    n(i+1)=n(i)+deltaT*(alpha_n(i)*(1-n(i))-beta_n(i)*n(i));
    m(i+1)=m(i)+deltaT*(alpha_m(i)*(1-m(i))-beta_m(i)*m(i));
    h(i+1)=h(i)+deltaT*(alpha_h(i)*(1-h(i))-beta_h(i)*h(i));
end



V=V-15;
plot(t,V);
ylabel('Voltage')
xlabel('time')
title('Voltage Over Time');

figure
p1=plot(t,gNa*(m.^3).*h);
title('Potassium Conductance')
figure
p2=plot(t,gK*(n.^4));
title('Sodium Conductance')
