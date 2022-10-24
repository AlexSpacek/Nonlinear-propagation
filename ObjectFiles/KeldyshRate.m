%           tic 
            n=1.76;
            m=1; %AU
            mSI=m*9.10938356*1e-31;
%             Eg=9/27.21; %ev to J 1.60218e-19
            EgSI=6.5*1.60218e-19;
            omega0SI=2*pi*2.997e8/(0.8e-6);
            E=linspace(1e9,3e10,1000);
            ySI=omega0SI*sqrt(m*9.10938356*1e-31*EgSI)./(1.60218e-19*E);
%             yAI=(omega0SI/4.13e16)*sqrt(m*9/27.21)./(E/5.1422e11);
            eSI=1.60218e-19;%*2.998e9
            omega0=omega0SI/4.13e16;
            h_barSI=1.054571800e-34; %SI
%             F2=2*1/2*abs(E/5.14220652e11).^2;
            I=n*2.65*1e-3*1/2*abs(E).^2;
            Icm=I*1e-4;
            Eg_barSI=EgSI+eSI^2*E.^2/(4*mSI*(omega0SI)^2);
            xSI=Eg_barSI/(h_barSI*omega0SI);
            PhiVar=sqrt(2*floor(xSI+1)-2*xSI);
            PhiVarSI=sqrt(2*floor(xSI+1)-2*xSI);
%             [~,indexMax]=max(E);
            intDawson = dawson(PhiVarSI);
%             intDawson=0.54;
            W_SI=(2*omega0SI/(9*pi))*(mSI*omega0SI/h_barSI)^(3/2)...
            .*intDawson.*exp(2*floor(xSI+1).*(1-(eSI^2*(E).^2)./(4*mSI*(omega0SI)^2*EgSI)))...
            .*(eSI^2*(E).^2./(16*mSI*EgSI*(omega0SI)^2)).^(floor(xSI+1)); 

            rho_nt=2.35e28; %fused silica m^-3
            mineSI=W_SI;
            mineCmSI=mineSI*1e-6;
            
%            toc 
            
            %Tunneling
% c = 3e8;
% Electron Charge [C]
% e = 1.6e-19;
% % Permittivity of Free Space [F/m]
% ep0 = 8.85e-12;
% % Planck Constant [J.s]
% hbar = 1.054*10^(-34);


% F = E;

Wtun1 = (2*EgSI)/(9*h_barSI*pi()^2) * ((mSI*EgSI)/h_barSI^2)^(3/2);
Wtun2 = ((eSI*h_barSI*E)./(mSI^(1/2)*EgSI^(3/2))).^(5/2);


Wtun3 = exp(-(pi()*mSI^(1/2)*EgSI^(3/2))./(2*eSI*h_barSI*E) .* (1-(mSI*omega0SI...
^2*EgSI)./(8*eSI^2*E.^2)));

Wtunnel = Wtun1 .* Wtun2 .*Wtun3; % [electrons/s/mˆ3]
%         toc   
figure
loglog(Icm,W_SI)
hold on
loglog(Icm,Wtunnel)        %Tunneling pro ySI=1.10;   
loglog(Icm,1e6*Icm.^7*4.3e-78/(7*h_barSI*omega0SI))

Wfull=KeldyshFullRate(omega0SI,mSI,EgSI,n,I);
loglog(Icm,Wfull)    
            
legend('MPI','TUN','Approx','Full Keldysh')            

omegaVec=linspace(2*pi*3e8/(600e-9),2*pi*3e8/(100e-9),1000);
absCoeff=real(sqrt(h_barSI*omegaVec-0.8*EgSI)./(h_barSI*omegaVec));
figure
plot(2*pi*3e8./omegaVec,absCoeff)
