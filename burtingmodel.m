%% SFO Model

% Burst Firing:
% For values of gNa=150 & gK=100 (in mS/cm^2) we get burst firing,
% gKS ranges from 2-5mS/cm^2 depending on the burst regime (B1 vs B2).

% Tonic Firing:
% Transition burst firing regime to tonic firing by either:
% 1. Increasing gK from 100mS/cm^2 to 280mS/cm^2.
% 2. Increasing gNa from from 150mS/cm^2 to 170mS/cm^2.

clear

%% Time:
dt = 0.01;                % Step Size for Eulers Method
timesec = 30;             % End Time (s)
tfinal = timesec*1000;    % End Time (ms)
Tstep =(0:dt:tfinal);     % Time Array (ms)

%% Initial Values:
V=-60;                    % Membrane Potential (mV)
nA=0;                     % IK Activation
mA=0.35;                  % IKS Activation
hA=0;                     % IKS Inactivation
mKA=0;                    % IA Activation
hKA=0;                    % IA Inactivation
sA=0;                     % INa Activation
kA=0;                     % INa Inactivation
mNaP=0;                   % INaP Activation
hNaP=0;                   % INaP Inactivation
mNCa=0;                   % ICa Activation

%% Cell Properties
% Cell Morphology:
CellSize = 10;                       % Model Cell Diameter in microns
CellArea = 4*pi*((CellSize/2).^2);   % Model Cell Area for Spherical Shape
% Cell Capacitance:
Cr = 5;                              % Real Cell Capacitance in pF
C = (Cr/CellArea)*100;               % Model Capacitance in uF/cm^2
% Cell Resistance:
Rr = 1;                              % Real Cell Resistance in gigaohms
R = (Rr * CellArea)*10;              % Model Cell Resistance in ohms*cm^2
% Reversal Potential (in mV):
Er = -65;                            % Reversal Potential for leak
ENa = 107;                           % Reversal Potential for Na+
EK = -88;                            % Reversal Potential for K+
ENSCC = -35;                         % Reversal Potential for NSCC
ECa = 120;                           % Reversal Potential for Ca2+
% Conductance (in mS/cm^2):
gLeak = (1/R)*1000;                  % Conductance of IL
gK = 100;                            % Conductance of IK
gKS = 3;                             % Conductance of IKS
gA = 3;                              % Conductance of IA
gNa = 150;                           % Conductance of INa
gNaP = 0.13;                         % Conductance of INaP
gNSCC = 0.2;                         % Conductance of INSCC
gCa = 0.3;                           % Conductance of ICa

%% Zeros Vectors
Vm=zeros(1,length(Tstep));
Iz=zeros(1,length(Tstep));
IN=zeros(1,length(Tstep));

%% Looping Code
% Solving differential equations with Eulers method:

tidx=0;
for z=Tstep
    tidx=tidx+1;

%% Injected Current
    % Where 1 pA = 0.3183 uA/cm^2 & 10 pA = 3.18 uA/cm^2
    I=0;

%% Time Constants (in ms):
    tau_nA = 7.2-(6.4/(1+exp((V+28.3)/-19.2)));  % Time Constant for IK activation
    tau_mA = 3000;                               % Time Constant for IKS activation
    tau_hA = 10;                                 % Time Constant for IKS inactivation
    tau_sA = 0.1;                                % Time Constant for INa activation
    tau_kA = 1;                                  % Time Constant for INa inactivation
    tau_mNaP = 5;                                % Time Constant for INaP activation
    tau_hNaP = 50;                               % Time Constant for INaP inactivation
    tau_mNCa = 5;                                % Time Constant for ICa activation
    tau_mKA = 5;                                 % Time Constant for IA activation
    tau_hKA = 30;                                % Time Constant for IA inactivation

%% Current Equations:
  % Potassium Currents:
    % Delayed-rectifier K+ Current (IK):
    nA0=1/(1+exp((V-2)/-8));
    nA = nA + dt*((nA0-nA)/tau_nA);
    IK = gK*nA^4*(V-EK);
    % Slow-activating K+ Current (IKS):
    mA0 = 1/(1+exp((V+44)/-18));
    mA = mA + dt*((mA0-mA)/tau_mA);
    hA0 = 1/(1+exp((V+61)/8));
    hA = hA + dt*((hA0-hA)/tau_hA);
    IKS = gKS*mA^3*hA*(V-EK);
    % Transient K+ Current (IA):
    mKA0 = 1/(1+exp((V+44)/-18));
    mKA = mKA + dt*((mKA0-mKA)/tau_mKA);
    hKA0 = 1/(1+exp((V+61)/8));
    hKA = hKA + dt*((hKA0-hKA)/tau_hKA);
    IA = gA*mKA^3*hKA*(V-EK);
  % Sodium Currents:
    % Transient Na+ Current (INa):
    sA0 = 1/(1+exp((V+31)/-6.11));
    sA = sA + dt*((sA0-sA)/tau_sA);
    kA0 = 1/(1+exp((V+62)/6.15));
    kA = kA + dt*((kA0-kA)/tau_kA);
    INa = gNa*sA^3*kA*(V-ENa);
    % Persisent Na+ Current (INaP):
    mNaP0 = 1/(1+exp((V+55)/-4));
    mNaP = mNaP + dt*((mNaP0-mNaP)/tau_mNaP);
    hNaP0 =1/(1+exp((V+45)/6.1));
    hNaP = hNaP + dt*((hNaP0-hNaP)/tau_hNaP);
    INaP = gNaP*mNaP^3*hNaP*(V-ENa);
  % Calcium Currents:
    % N-Type Ca2+ Current (ICa):
    mNCa0 = 1/(1+exp((V+14)/-5.8));
    mNCa = mNCa + dt*((mNCa0-mNCa)/tau_mNCa);
    ICa = gCa*mNCa^2*(V-ECa);
  % Leak Current (IL):
    IL = gLeak*(V-Er);
  % Non-selective Cation Current (INSCC):
    INSCC = gNSCC*(V-ENSCC);
  % Noise (IN):
    stdnoise = 9;          % Standard Deviation for Noise
    IN=stdnoise*randn;     % Noise Current (Gaussian Distribution)

%% Voltage Calculation:
    V = V + ((dt/C)*(I + IN - IL - INa - IK - INSCC - INaP - IA - ICa - IKS));

%% Voltage Array:
    Vm(tidx)= V;

end

%% Membrane Potential vs Time Plot:
figure('Renderer', 'painters', 'Position', [100 100 1200 500])
sh(1)=subplot(1,1,1);
plot(Tstep,Vm,'k-')
% Figure Settings:
box off
ax = gca;
ax.FontSize = 14; 
xlabel('Time (s)','FontSize',14)
ylabel('Membrane Potential (mV)','FontSize',14)
xt=arrayfun(@num2str,get(gca,'xtick')/1000,'un',0);
set(gca,'xticklabel',xt)