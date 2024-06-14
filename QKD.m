clc; clear; close all;
% Data SPAD
tBob=10^(-3/10); % Bob Internal losses
KmTot = 500;
KmVect = 0:0.1:KmTot;
tAB= 10.^((-0.2/10).*KmVect);
etaD=0.2 % detector efficiency 
eta= tBob*etaD.*tAB % overall detection efficiency
pdc= 10^(-5); %Dark counts
emis = 10^(-2) %Misalignment for which a photon hits the erroneous detector

% SPAD param.
Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Q1=Y1; % Gain (No Decoy State)
QBER = (pdc+emis.*eta)./(eta+(1-eta).*Y0);
% Superconductor data
tBobs= 10^(-3/10);
etaDs = 0.9;
etas = tBobs*etaDs*tAB;
pdcs = 10^(-9); 
emiss = 10^(-2);

% Superconductor param.
Y0s = 2*pdcs-pdcs^2;
Y1s = etas;
Q1s=Y1s;
QBERs = emis+(pdcs./etas);
% SKR and Plotting
% For QBERs the fact we are using other formula implies we need to constrain the
% limit for QBERs in order to not surpass the Shannon Entropy H(QBERs)<=1, so => QBERs<= 0.5

ZrosQBERs = find(QBERs>=0.5);
QBERs(1,ZrosQBERs)=0.5;

ZrosQBER = find(QBER>=0.5);
QBER(1,ZrosQBER)=0.5;

SKR = 0.5.*Q1.*(1-2.*H(QBER)); %SKR for SPAD (Single Photon Avalanche Detection)
SKRs = 0.5.*Q1s.*(1-2.*H(QBERs));% SKR for Superconductor

ZrosSKRs = find(SKRs<=0);% Fixing negative values for high values of KmVect
SKRs(1,ZrosSKRs)=0; % Swap of these negative values for 0
ZrosSKR = find(SKR<=0);
SKR(1,ZrosSKR)=0;
% Visual Plot Refinement 
GraphSPAD= find(SKR==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKR(1,GraphSPAD)=1e-12;

GraphSuperconductor= find(SKRs==0);
SKRs(1,GraphSuperconductor)=1e-12;
for t = 51:length(SKR)
    if SKR(t) <= 0.5*SKR(t-50);
        distance_fall = KmVect(t);
        valueSKR = SKR(t);
        break;
    end
end
for t = 51:length(SKRs)
    if SKRs(t) <= 0.5.*SKRs(t-50)
        distance_fallS = KmVect(t);
        valueSKRs = SKRs(t);
        break;
    end
end
threedBfallSKR  = valueSKR.*ones(1,5001)
threedBfallSKRs = valueSKRs.*ones(1,5001)

% Plotting SKR VS Km
semilogy(KmVect,SKRs)% Log representation of the y axis
hold on
semilogy(KmVect,SKR)
axis([0 500 10e-11 1])
lgd= legend('SSPD', 'SPAD');
fontsize(lgd,9,'points')
xlabel("Fiber Length [Km]")
ylabel("Key rate per pulse [bit/pulse]")

semilogy(KmVect,threedBfallSKR,"--r")
semilogy(KmVect,threedBfallSKRs,"--r")
semilogy(distance_fall,valueSKR,'rx','MarkerSize',10)
semilogy(distance_fallS,valueSKRs,'rx','MarkerSize',10)
lgd= legend('SSPD', 'SPAD');
fontsize(lgd,9,'points')
xlabel("Fiber Length [Km]")
ylabel("Key rate per pulse [bit/pulse]")

hold off
etaDper = 0.000001:0.00001:1
tABper= 1 %0.00003162277 we begin the simulations with 45dB losses, my computer cannot without adding this losses, so much data.
 etaper= 1.*etaDper.*tABper;
 figure();
 pdcper= 10^(-5); %Dark counts SPAD
 emisper = 10^(-2); %Misalignment (0%,1%,...,5%)
 Y0per= 2*pdcper-pdcper^2; %Yield void state
 Y1per = Y0per+etaper-Y0per.*etaper; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
 Q1per=Y1per; % Gain (No Decoy State)
 QBERper = (pdcper+emisper.*etaper)./(etaper+(1-etaper).*Y0per);
 ZrosQBERper = find(QBERper>=0.5);
 QBERper(1,ZrosQBERper)=0.5;
 SKRper = 0.5.*Q1per.*(1-2.*H(QBERper));
 ZrosSKRper = find(SKRper<=0);
 SKRper(1,ZrosSKRper)=0;
 GraphSPADper= find(SKRper==0); % This makes more visual the result, even though its not "real" but we can approximate it
 SKRper(1,GraphSPADper)=1e-12;
            dB = 10.*log10(etaper); % From linear to dB
            semilogy(-dB,SKRper)
            axis([0 100 10e-12 1])
            xlabel("Overall attenuation (-η)[dB]")
            ylabel("Key rate per pulse [bit/pulse]")
%% (SPAD)Detection efficency VS SKR, param. emis(misalignment), dark count and perfect Bob detector with no losses.

%Data, param. --> misalingment(subplots) and dark counts (plots)

tBobper=10^(-0/10); % 
tABper= 10.^((-0.2/10).*0); % We do the simulations with no attenuation so km=0;
etaDper = 0.00000001:0.000001:1; % detector efficiency
% We devide the subplots in to loops because we want two figures with the
% subplots
figure();
  for j = 2:1:5 % Dark count (10^(-2),10^(-3),...,10^(-5))
      order_index = j-1; % aux. variable
      subplot(2, 2, order_index) ;
        for i = 0:5 % Misalignment (0%,1%,...,5%)
            etaper= tBobper.*etaDper.*tABper;
            pdcper= 10^(-j); %Dark counts SPAD
            emisper = i*10^(-2); %Misalignment (0%,1%,...,5%)
        
        
            Y0per= 2*pdcper-pdcper^2; %Yield void state
            Y1per = Y0per+etaper-Y0per.*etaper; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice)
            Q1per=Y1per; % Gain (No Decoy State)
            QBERper = (pdcper+emisper.*etaper)./(etaper+(1-etaper).*Y0per);
            ZrosQBERper = find(QBERper>=0.5);
            QBERper(1,ZrosQBERper)=0.5;
            SKRper = 0.5.*Q1per.*(1-2.*H(QBERper));
            ZrosSKRper = find(SKRper<=0);
            SKRper(1,ZrosSKRper)=0;
            GraphSPADper= find(SKRper==0); % This makes more visual the result, even though its not "real" but we can approximate it
            SKRper(1,GraphSPADper)=1e-12;
            dB = -10.*log10(etaper); % From linear to dB
            semilogy(dB,SKRper,'DisplayName',"Misalignment" + i +"%")
            axis([15 45 10e-8 1])
             
            hold on
            
        end
      hold off
      
      
     lgd= legend('Location','northeast');
     fontsize(lgd,7,'points')
      odd_check= mod(order_index,2); % Odd number to avoid cluttering the graph with unnecessary text
      if odd_check == 1
        xlabel("Overall attenuation (-η)[dB]")
        ylabel("Key rate per pulse [bit/pulse]")
      else
          xlabel("Overall attenuation (-η)[dB]")
      end
      title(sprintf('Dark Count order (-%d)',j))     
      sgtitle('SPAD - SKR vs Efficiency: Dark Count Rate Impact')
     
  end
  %%
  figure();
   for j = 6:1:9 % Dark count (10^(-6),10^(-7),...,10^(-9))
      order_index = j-5; % To represent...
      subplot(2, 2, order_index) ;
        for i = 0:5 % misalignment (0%,1%,...,5%)
            tABper= 0.00003162277 % we begin the simulations with 45dB losses, my computer cannot without adding this losses, so much data.
            etaper= tBobper.*etaDper.*tABper;
            pdcper= 10^(-j); %Dark counts SPAD
            emisper = i*10^(-2); 
        
            % SPAD CASE
        
            Y0per= 2*pdcper-pdcper^2; %Yield void state
            Y1per = Y0per+etaper-Y0per.*etaper; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice)
            Q1per=Y1per; % Gain (No Decoy State)
            QBERper = (pdcper+emisper.*etaper)./(etaper+(1-etaper).*Y0per);
            ZrosQBERper = find(QBERper>=0.5);
            QBERper(1,ZrosQBERper)=0.5;
            SKRper = 0.5.*Q1per.*(1-2.*H(QBERper));
            ZrosSKRper = find(SKRper<=0);
            SKRper(1,ZrosSKRper)=0;
            GraphSPADper= find(SKRper==0); % This makes more visual the result, even though its not "real" but we can approximate it
            SKRper(1,GraphSPADper)=1e-12;
            dB = 10.*log10(etaper); % From linear to dB
            semilogy(-dB,SKRper,'DisplayName',"Misalignment" + i +"%")
            axis([55 85 10e-11 10e-4])
             
            hold on
            
        end
      hold off
      
      
      lgd= legend('Location','northeast');
      fontsize(lgd,7,'points')
      odd_check= mod(order_index,2);
      if odd_check == 1
        xlabel("Overall attenuation (-η)[dB]")
        ylabel("Key rate per pulse [bit/pulse]")
      else
          xlabel("Overall attenuation (-η)[dB]")
      end
      title(sprintf('Dark Count order (-%d)',j))     
      sgtitle('SKR vs Efficiency: Dark Count Rate Impact')
     
   end

