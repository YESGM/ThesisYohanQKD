%% Decoy state case
clc; clear;
% Data SPAD
tBob=10^(-3/10); % Bob Internal losses
KmTot = 500;
KmVect = 0:0.1:KmTot;
tAB= 10.^((-0.2/10).*KmVect);
etaD=0.2 % detector efficiency 
eta= tBob*etaD.*tAB % overall detection efficiency
pdc= 10^(-5); %Dark counts
emis = 10^(-2) %Misalignment for which a photon hits the erroneous detector
mu=0.3; % expected photon number by Alice(original signal)
nu=0.05; % (decoy state v)
emis=10^(-2);

%-----------SKR for attenuated lasers with decoy state (mu,nu(random),0)------------
%
% SPAD param.
Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Qmu=Y0+1-exp(-eta.*mu);
Qnu=Y0+1-exp(-eta.*nu);
QBERmu = (0.5*Y0+emis.*(1-exp(-eta.*mu)))./Qmu;
QBERnu = (0.5*Y0+emis.*(1-exp(-eta.*nu)))./Qnu;
factor_low_bound= (Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
YLv0=(mu./(mu.*nu-nu.^2)).*factor_low_bound;
QLv0=((mu.^2.*exp(-mu))./(mu.*nu-nu.^2)).*(Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
e1Lv0 = (QBERnu.*Qnu.*exp(nu)-0.5.*Y0)./(YLv0.*nu);

SKRlv0 = 0.5.*(-Qmu.*H(QBERmu)+QLv0.*(1-H(e1Lv0))); %SKR for SPAD (Single Photon Avalanche Detection)


ZrosSKRlv0 = find(SKRlv0<=0);

SKRlv0(1,ZrosSKRlv0)=0;

GraphSPADlv0= find(SKRlv0==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKRlv0(1,GraphSPADlv0)=1e-12;

figure();
semilogy(KmVect,SKRlv0)
hold on

%-----------SKR for attenuated lasers with decoy state (mu,nu(optimized),0)------------
clear mu
syms mu 
t=solve(((H(emis)./(1-H(emis))))==exp(-mu).*(1-mu),mu); % Optimization of mu
comp=((H(emis)./(1-H(emis))));
mu= double(t); % Optimal value of mu

Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Qmu=Y0+1-exp(-eta.*mu);
Qnu=Y0+1-exp(-eta.*nu);
Q1=Y1.*mu.*exp(-mu);
e1=(0.5*Y0+emis.*eta)./Y1;
QBERmu = (0.5*Y0+emis.*(1-exp(-eta.*mu)))./Qmu;
QBERnu = (0.5*Y0+emis.*(1-exp(-eta.*nu)))./Qnu;
factor_low_bound= (Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
YLv0=(mu./(mu.*nu-nu.^2)).*factor_low_bound;
QLv0=((mu.^2.*exp(-mu))./(mu.*nu-nu.^2)).*(Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
e1Lv0 = (QBERnu.*Qnu.*exp(nu)-0.5.*Y0)./(YLv0.*nu);

SKRlv0 = 0.5.*(-Qmu.*H(QBERmu)+QLv0.*(1-H(e1Lv0))); %SKR for SPAD (Single Photon Avalanche Detection) using gain low bound and upper bound for the quantum bit error rate, using mu optimized.


ZrosSKRlv0 = find(SKRlv0<=0);

SKRlv0(1,ZrosSKRlv0)=0;

GraphSPADlv0= find(SKRlv0==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKRlv0(1,GraphSPADlv0)=1e-12;
semilogy(KmVect,SKRlv0, '--')

%-----------SKR for attenuated lasers with non-decoy state-----------------

clear mu
mu=0.4; %random mu

Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Qmu=Y0+1-exp(-eta.*mu);
Qnu=Y0+1-exp(-eta.*nu);
Q1=Y1.*mu.*exp(-mu);
e1=(0.5*Y0+emis.*eta)./Y1;
QBERmu = (0.5*Y0+emis.*(1-exp(-eta.*mu)))./Qmu;
QBERnu = (0.5*Y0+emis.*(1-exp(-eta.*nu)))./Qnu;
factor_low_bound= (Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
YLv0=(mu./(mu.*nu-nu.^2)).*factor_low_bound;
QLv0=((mu.^2.*exp(-mu))./(mu.*nu-nu.^2)).*(Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
e1Lv0 = (QBERnu.*Qnu.*exp(nu)-0.5.*Y0)./(YLv0.*nu);

SKRnodec=0.5.*(-Qmu.*H(QBERmu)+Q1.*(1-H(e1))) %Multiphoton SKR without decoy state, omega =/ 1.

ZrosSKRnodec = find(SKRnodec<=0);

SKRnodec(1,ZrosSKRnodec)=0;

GraphSPADnodec= find(SKRnodec==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKRnodec(1,GraphSPADnodec)=1e-12;

semilogy(KmVect,SKRnodec)

%-----------SKR for the initial single-photon sources----------------------

Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Q1=Y1; % Gain (No Decoy State)
QBER = (pdc+emis.*eta)./(eta+(1-eta).*Y0);

ZrosQBER = find(QBER>=0.5);
QBER(1,ZrosQBER)=0.5;

SKR = 0.5.*Q1.*(1-2.*H(QBER)); % SKR with ideal single-photon sources
ZrosSKR = find(SKR<=0);
SKR(1,ZrosSKR)=0;
% Visual Plot Refinement 
GraphSPAD= find(SKR==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKR(1,GraphSPAD)=1e-12;

semilogy(KmVect,SKR)
axis([0 500 10e-12 1])
legend
xlabel("Fiber Length [Km]")
ylabel("Key rate per pulse [bit/pulse]")
legend('Decoy State random mu','Decoy State optimal mu','Non Decoy State', 'Single-photon source')
hold off
 lgd= legend('Location','northeast');
  fontsize(lgd,9,'points')
title('SKR VS Distance - SPAD') 

  hold off

%--------------------SUPERCONDUCTOR CASE-----------------------------------

% Decoy state case
clc; clear;
% Data Superconductor
tBob=10^(-3/10); % Bob Internal losses
KmTot = 500;
KmVect = 0:0.1:KmTot;
tAB= 10.^((-0.2/10).*KmVect);
etaD=0.9; % detector efficiency 
eta= tBob*etaD.*tAB % overall detection efficiency
pdc= 10^(-9); %Dark counts
emis = 10^(-2); %Misalignment for which a photon hits the erroneous detector
mu=0.3; % expected photon number by Alice(original signal)
nu=0.05; % (decoy state v)
emis=10^(-2);

%-----------SKR for attenuated lasers with decoy state (mu,nu(random),0)------------

% Superconductor param.
Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Qmu=Y0+1-exp(-eta.*mu);
Qnu=Y0+1-exp(-eta.*nu);
QBERmu = (0.5*Y0+emis.*(1-exp(-eta.*mu)))./Qmu;
QBERnu = (0.5*Y0+emis.*(1-exp(-eta.*nu)))./Qnu;
factor_low_bound= (Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
YLv0=(mu./(mu.*nu-nu.^2)).*factor_low_bound;
QLv0=((mu.^2.*exp(-mu))./(mu.*nu-nu.^2)).*(Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
e1Lv0 = (QBERnu.*Qnu.*exp(nu)-0.5.*Y0)./(YLv0.*nu);

SKRlv0 = 0.5.*(-Qmu.*H(QBERmu)+QLv0.*(1-H(e1Lv0))); %SKR for attenuated weak coherent pulses, using decoy state random mu


ZrosSKRlv0 = find(SKRlv0<=0);

SKRlv0(1,ZrosSKRlv0)=0;

GraphSPADlv0= find(SKRlv0==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKRlv0(1,GraphSPADlv0)=1e-12;

figure();
semilogy(KmVect,SKRlv0)
hold on

%-----------SKR for attenuated lasers with decoy state (mu,nu(optimized),0)------------

clear mu
syms mu 
t=solve(((H(emis)./(1-H(emis))))==exp(-mu).*(1-mu),mu); % Optimization of mu
comp=((H(emis)./(1-H(emis))));
mu= double(t);

Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Qmu=Y0+1-exp(-eta.*mu);
Qnu=Y0+1-exp(-eta.*nu);
Q1=Y1.*mu.*exp(-mu);
e1=(0.5*Y0+emis.*eta)./Y1;
QBERmu = (0.5*Y0+emis.*(1-exp(-eta.*mu)))./Qmu;
QBERnu = (0.5*Y0+emis.*(1-exp(-eta.*nu)))./Qnu;
factor_low_bound= (Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
YLv0=(mu./(mu.*nu-nu.^2)).*factor_low_bound;
QLv0=((mu.^2.*exp(-mu))./(mu.*nu-nu.^2)).*(Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
e1Lv0 = (QBERnu.*Qnu.*exp(nu)-0.5.*Y0)./(YLv0.*nu);

SKRlv0 = 0.5.*(-Qmu.*H(QBERmu)+QLv0.*(1-H(e1Lv0))); %SKR for Superconductor using gain low bound and upper bound for the quantum bit error rate, optimized mu


ZrosSKRlv0 = find(SKRlv0<=0);

SKRlv0(1,ZrosSKRlv0)=0;

GraphSPADlv0= find(SKRlv0==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKRlv0(1,GraphSPADlv0)=1e-12;
semilogy(KmVect,SKRlv0, '--')
axis([0 500 10e-12 1])
legend
xlabel("Fiber Length [Km]")
ylabel("Key rate per pulse [bit/pulse]")
legend('decoy state with random mu','decoy state optimal mu')

%-----------SKR for attenuated lasers with non-decoy state------------------

clear mu
mu=0.4; %random mu

Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Qmu=Y0+1-exp(-eta.*mu);
Qnu=Y0+1-exp(-eta.*nu);
Q1=Y1.*mu.*exp(-mu);
e1=(0.5*Y0+emis.*eta)./Y1;
QBERmu = (0.5*Y0+emis.*(1-exp(-eta.*mu)))./Qmu;
QBERnu = (0.5*Y0+emis.*(1-exp(-eta.*nu)))./Qnu;
factor_low_bound= (Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
YLv0=(mu./(mu.*nu-nu.^2)).*factor_low_bound;
QLv0=((mu.^2.*exp(-mu))./(mu.*nu-nu.^2)).*(Qnu.*exp(nu)-Qmu.*exp(mu).*((nu.^2)./mu.^2)-(((mu.^2 - nu.^2)./mu.^2).*Y0));
e1Lv0 = (QBERnu.*Qnu.*exp(nu)-0.5.*Y0)./(YLv0.*nu);

SKRnodec=0.5.*(-Qmu.*H(QBERmu)+Q1.*(1-H(e1)));

ZrosSKRnodec = find(SKRnodec<=0);

SKRnodec(1,ZrosSKRnodec)=0;

GraphSPADnodec= find(SKRnodec==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKRnodec(1,GraphSPADnodec)=1e-12;

semilogy(KmVect,SKRnodec)

%-----------SKR for the initial single-photon sources----------------------------------

Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Q1=Y1; % Gain (No Decoy State)
QBER = (pdc+emis.*eta)./(eta+(1-eta).*Y0);

ZrosQBER = find(QBER>=0.5);
QBER(1,ZrosQBER)=0.5;

SKR = 0.5.*Q1.*(1-2.*H(QBER)); % SKR with ideal single-photon sources
ZrosSKR = find(SKR<=0);
SKR(1,ZrosSKR)=0;
% Visual Plot Refinement 
GraphSPAD= find(SKR==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKR(1,GraphSPAD)=1e-12;

semilogy(KmVect,SKR)
axis([0 500 10e-12 1])
legend
xlabel("Fiber Length [Km]")
ylabel("Key rate per pulse [bit/pulse]")
legend('Decoy State random mu','Decoy State optimal mu','Non Decoy State', 'Single-photon source', 'Single-photon source')
 lgd= legend('Location','northeast');
  fontsize(lgd,9,'points')

title('SKR VS Distance - Superconductor') 
hold off