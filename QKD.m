clc; clear; close all;
%% Data SPAD
tBob=10^(-3/10); % Bob Internal losses
KmTot = 500;
KmVect = 0:0.1:KmTot;
tAB= 10.^((-0.2/10).*KmVect);
etaD=0.2; % detector efficiency 
eta= tBob*etaD.*tAB % overall detection efficiency
pdc= 10^(-5); %Dark counts
emis = 10^(-2) %Misalignment for which a photon hits the erroneous detector

% SPAD param.
Y0= 2*pdc-pdc^2; %Yield void state
Y1 = Y0+eta-Y0.*eta; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
Q1=Y1; % Gain (No Decoy State)
QBER = (pdc+emis.*eta)./(eta+(1-eta).*Y0);
%% Superconductor data
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
%% SKR and Plotting
% For QBERs the fact we are using other formula implies we need to constrain the
% limit for QBERs in order to not surpass the Shannon Entropy H(QBERs)<=1, so => QBERs<= 0.5

ZrosQBERs = find(QBERs>=0.5);
QBERs(1,ZrosQBERs)=0.5;

ZrosQBER = find(QBER>=0.5);
QBER(1,ZrosQBER)=0.5;

SKR = 0.5.*Q1.*(1-2*H(QBER)); %SKR for SPAD (Single Photon Avalanche Detection)
SKRs = 0.5.*Q1s.*(1-2*H(QBERs));% SKR for Superconductor

ZrosSKRs = find(SKRs<=0);% Fixing negative values for high values of KmVect
ZrosSKR = find(SKR<=0);
SKRs(1,ZrosSKRs)=0; % Swap of these negative values for 0
SKR(1,ZrosSKR)=0;
%% Visual Plot Refinement 
GraphSPAD= find(SKR==0); % This makes more visual the result, even though its not "real" but we can approximate it
SKR(1,GraphSPAD)=1e-12;

GraphSuperconductor= find(SKRs==0);
SKRs(1,GraphSuperconductor)=1e-12;
%% Plotting SKR VS Km
semilogy(KmVect,SKRs)% Log representation of the y axis
hold on
semilogy(KmVect,SKR)
axis([0 400 10e-12 1])
legend
xlabel("Fiber Length [Km]")
ylabel("Key rate per pulse [bit/pulse]")
legend('Superconductor','SPAD')
hold off
%% (SPAD)Detection efficency VS SKR, param. emis(misalignment), dark count and perfect Bob detector with no losses.

%Data, param. --> misalingment(subplots) and dark counts (plots)

tBobper=10^(-0/10); % 
% KmTotper = 100;
% KmVectper = 0:1:KmTotper;
tABper= 10.^((-0.2/10).*0); % We do the simulations with no attenuation so km=0;
etaDper = 0:0.01:1; % detector efficiency
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
        
            % SPAD CASE
        
            Y0per= 2*pdcper-pdcper.^2; %Yield void state
            Y1per = Y0per+etaper-Y0per.*etaper; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
            Q1=Y1per; % Gain (No Decoy State)
            QBERper = (pdcper+emisper.*etaper)./(etaper+(1-etaper).*Y0per);
        
            ZrosQBERper = find(QBERper>=0.5);
            QBERper(1,ZrosQBERper)=0.5;
        
            SKRper = 0.5.*Q1.*(1-2*H(QBERper));
            ZrosSKRper = find(SKRper<=0);
            SKRper(1,ZrosSKRper)=0;
            dB = mag2db(abs(etaper)); % From linear to dB
                if i==0
                    plot(dB,SKRper,'DisplayName',"Misalignment " + i + "%")
                else
                    plot(dB,SKRper,'DisplayName',"Misalignment " + i + "0%")
                end
            hold on
            
        end
      hold off
      
    
      lgd=legend('Location','northwest');
      fontsize(lgd,7,'points')
      odd_check= mod(order_index,2); % Odd number to avoid cluttering the graph with unnecessary text
      if odd_check == 1
        xlabel("η (Total detection efficiency)[dB]")
        ylabel("Key rate per pulse [bit/pulse]")
      else
          xlabel("η (Total detection efficiency)[dB]")
      end
      title(sprintf('Dark Count order (-%d)',j))     
      sgtitle('SPAD - SKR vs Efficiency: Dark Count Rate Impact')
     
  end
  figure();
   for j = 6:1:9 % Dark count (10^(-6),10^(-7),...,10^(-9))
      order_index = j-5; % To represent...
      subplot(2, 2, order_index) ;
        for i = 0:5 % misalignment (0%,1%,...,5%)

            etaper= tBobper.*etaDper.*tABper;
            pdcper= 10^(-j); %Dark counts SPAD
            emisper = i*10^(-2); 
        
            % SPAD CASE
        
            Y0per= 2*pdcper-pdcper.^2; %Yield void state
            Y1per = Y0per+etaper-Y0per.*etaper; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice) 
            Q1=Y1per; % Gain (No Decoy State)
            QBERper = (pdcper+emisper.*etaper)./(etaper+(1-etaper).*Y0per);
        
            ZrosQBERper = find(QBERper>=0.5);
            QBERper(1,ZrosQBERper)=0.5;
        
            SKRper = 0.5.*Q1.*(1-2*H(QBERper));
            ZrosSKRper = find(SKRper<=0);
            SKRper(1,ZrosSKRper)=0;
            dB = mag2db(abs(etaper));
                if i==0
                    plot(dB,SKRper,'DisplayName',"Misalignment " + i + "%")
                else
                    plot(dB,SKRper,'DisplayName',"Misalignment " + i + "0%")
                end
            hold on
            
        end
      hold off
      
    

      lgd=legend('Location','northwest');
      fontsize(lgd,7,'points')
      odd_check= mod(order_index,2);
      if odd_check == 1
        xlabel("η (Total detection efficiency)[dB]")
        ylabel("Key rate per pulse [bit/pulse]")
      else
          xlabel("η (Total detection efficiency)[dB]")
      end
      title(sprintf('Dark Count order (-%d)',j))     
      sgtitle('SPAD - SKR vs Efficiency: Dark Count Rate Impact')
     
   end
%% (Superconductor)Detection efficency VS SKR, param. emis(misalignment), dark count and perfect Bob detector with no losses.

%Data, param. --> misalingment(subplots) and dark counts (plots)
clear 
clc

tBobper=10^(-0/10); % 
% KmTotper = 100;
% KmVectper = 0:1:KmTotper;
tABper= 10.^((-0.2/10).*0); % We do the simulations with no attenuation so km=0;
etaDper = 0:0.01:1; % detector efficiency
% We devide the subplots in to loops because we want two figures with the
% subplots
figure()
  for j = 2:1:5 % Dark count (10^(-2),10^(-3),...,10^(-5))
      order_index = j-1; % aux. variable
      subplot(2, 2, order_index) ;
        for i = 0:5 % Misalignment (0%,1%,...,5%)

            etaper= tBobper.*etaDper.*tABper;
            pdcper= 10^(-j); %Dark counts Superconductor
            emisper = i*10^(-2); %Misalignment (0%,1%,...,5%)
        
            % Superconductor CASE
        
            Y0per = 2*pdcper-pdcper^2;
            Y1per = etaper;
            Q1per=Y1per;
            QBERper = emisper+(pdcper./etaper);

            ZrosQBERper = find(QBERper>=0.5);
            QBERper(1,ZrosQBERper)=0.5;
        
            SKRper = 0.5.*Q1per.*(1-2*H(QBERper));
            ZrosSKRper = find(SKRper<=0);
            SKRper(1,ZrosSKRper)=0;
            dB = mag2db(abs(etaper)); % From linear to dB
                if i==0
                    plot(dB,SKRper,'DisplayName',"Misalignment " + i + "%")
                else
                    plot(dB,SKRper,'DisplayName',"Misalignment " + i + "0%")
                end
            hold on
            
        end
      hold off
      
    
      lgd=legend('Location','northwest');
      fontsize(lgd,7,'points')
      odd_check= mod(order_index,2); % Odd number to avoid cluttering the graph with unnecessary text
      if odd_check == 1
        xlabel("η (Total detection efficiency)[dB]")
        ylabel("Key rate per pulse [bit/pulse]")
      else
          xlabel("η (Total detection efficiency)[dB]")
      end
      title(sprintf('Dark Count order (-%d)',j))     
      sgtitle('Superconductor - SKR vs Efficiency: Dark Count Rate Impact')
     
  end
  figure();
   for j = 6:1:9 % Dark count (10^(-6),10^(-7),...,10^(-9))
      order_index = j-5; % To represent...
      subplot(2, 2, order_index) ;
        for i = 0:5 % misalignment (0%,1%,...,5%)

            etaper= tBobper.*etaDper.*tABper;
            pdcper= 10^(-j); %Dark counts SPAD
            emisper = i*10^(-2); 
        
            % Superconductor CASE
        
            Y0per = 2*pdcper-pdcper^2;
            Y1per = etaper;
            Q1per=Y1per;
            QBERper = emisper+(pdcper./etaper);

            ZrosQBERper = find(QBERper>=0.5);
            QBERper(1,ZrosQBERper)=0.5;
        
            SKRper = 0.5.*Q1per.*(1-2*H(QBERper));
            ZrosSKRper = find(SKRper<=0);
            SKRper(1,ZrosSKRper)=0;
            dB = mag2db(abs(etaper));
                if i==0
                    plot(dB,SKRper,'DisplayName',"Misalignment " + i + "%")
                else
                    plot(dB,SKRper,'DisplayName',"Misalignment " + i + "0%")
                end
            hold on
            
        end
      hold off
      
    

      lgd=legend('Location','northwest');
      fontsize(lgd,7,'points')
      odd_check= mod(order_index,2);
      if odd_check == 1
        xlabel("η (Total detection efficiency)[dB]")
        ylabel("Key rate per pulse [bit/pulse]")
      else
          xlabel("η (Total detection efficiency)[dB]")
      end
      title(sprintf('Dark Count order (-%d)',j))     
      sgtitle('Superconductor - SKR vs Efficiency: Dark Count Rate Impact')
     
   end

 %% Comparative SKR VS Efficiency SPAD VS Superconductor (Dark Count = 10^-5 and 10^-9, Misalignment = 10^-2 for both cases) with perfect BOB detector
 clear
 clc
   tBobper=10^(-0/10); % 
   % KmTotper = 100;
   % KmVectper = 0:1:KmTotper;
   tABper= 10.^((-0.2/10).*0); % We do the simulations with no attenuation so km=0;
   etaDper = 0:0.01:1; % detector efficiency
   etaper= tBobper.*etaDper.*tABper;
   pdcper= 10^(-5); %Dark counts SPAD
   emisper = 10^(-2); 
        
   % SPAD CASE

   Y0per= 2*pdcper-pdcper.^2; %Yield void state
   Y1per = Y0per+etaper-Y0per.*etaper; %Yield ( Bob detects a photon given a n-photon signal is emitted by Alice)
   Q1=Y1per; % Gain (No Decoy State)
   QBERper = (pdcper+emisper.*etaper)./(etaper+(1-etaper).*Y0per);

   ZrosQBERper = find(QBERper>=0.5);
   QBERper(1,ZrosQBERper)=0.5;

   SKRper = 0.5.*Q1.*(1-2*H(QBERper));
   ZrosSKRper = find(SKRper<=0);
   SKRper(1,ZrosSKRper)=0;
   dB = mag2db(abs(etaper));

   %Plotting
    plot(dB,SKRper,'DisplayName','SPAD')
    lgd=legend('Location','northwest');
    fontsize(lgd,12,'points')
    xlabel("η (Total detection efficiency)[dB]")
    ylabel("Key rate per pulse [bit/pulse]")
    hold on
clear; 
clc;
   tBobper=10^(-0/10); % Bob no internal losses
   % KmTotper = 100;
   % KmVectper = 0:1:KmTotper;
   tABper= 10.^((-0.2/10).*0); % We do the simulations with no attenuation so km=0;
   etaDper = 0:0.01:1; % detector efficiency
   etaper= tBobper.*etaDper.*tABper;
   pdcper= 10^(-9); %Dark counts Superconductor
   emisper = 10^(-2); 
        
   % Superconductor CASE

   Y0per = 2*pdcper-pdcper^2;
   Y1per = etaper;
   Q1per=Y1per;
   QBERper = emisper+(pdcper./etaper);

   ZrosQBERper = find(QBERper>=0.5);
   QBERper(1,ZrosQBERper)=0.5;

   SKRper = 0.5.*Q1per.*(1-2*H(QBERper));
   ZrosSKRper = find(SKRper<=0);
   SKRper(1,ZrosSKRper)=0;
   dB = mag2db(abs(etaper));

   %Plotting
    plot(dB,SKRper,'DisplayName','Superconductor')
    lgd=legend('Location','northwest');
    fontsize(lgd,12,'points')
    xlabel("η (Total detection efficiency)[dB]")
    ylabel("Key rate per pulse [bit/pulse]")
    hold off
