classdef Field < handle
    properties
        Data
        CentralWavelength
        FWHM
        Width
        WidthVector=zeros(1,1000)
        WidthVector2=zeros(1,1000)
        WidthVector3=zeros(1,1000)
        InitialEnergy
        Energy 
        MinWavelength
        maxSpace
        Size
        DispersionCoefficients 
        TimeVector
        FrequencyVector
        WavelengthVector
        SpaceVector
        SpaceFrequencyVector
        MeshSpace
        MeshTime
        WaveNumber
        CentralFrequency
        HenkelInputTransform
        HenkelOuputTransform
        HenkelC
        InitialSpectrum
        ZAxis
        bIntegral
        TimePropagate    
        TimeEnvelopePropagate
        SpectrumPropagate
        AverageSpectrumPropagate
        AverageSpectrumPropagate2
        RhoPropagate
        PrecalculatedE
        PrecalculatedW
        IndexSpectralCutoff
        
    end
    properties (Dependent)
        TemporalEnvelope
        SpaceEnvelope
        FrequencySpectrum
        SpaceFrequencySpectrum
    end    
    methods
        function thisField=Field(InitialParameters)              
            dimension=InitialParameters.Size;
            centralWavelength=InitialParameters.CentralWavelength;
            minWavelength=InitialParameters.MinWavelength;
            maxSpace=InitialParameters.MaxSpace;
            FWHM=InitialParameters.FWHM;
            width=InitialParameters.Width;
            energy=InitialParameters.Energy;
            thisField.InitialEnergy=energy;
            dispersionCoefficients=InitialParameters.DispersionCoefficients;
            thisField.ZAxis=zeros(1,1000);
            if length(fieldnames(InitialParameters)) == 8
                if isequal(size(dimension),[1 2])
                    thisField.Size=dimension;
                else
                    error('Size of the field matrix must be given for instance as [101 2001]')
                end    
                if isnumeric(centralWavelength)
                    thisField.CentralWavelength=centralWavelength;
                    thisField.WaveNumber=2*pi/centralWavelength;
                    thisField.CentralFrequency=2.9979e8/centralWavelength;
                else
                    error('Wavelength must be numeric')
                end    
                if isnumeric(minWavelength)
                    thisField.MinWavelength=minWavelength;

                else
                    error('Minimal wavelength must be numeric')
                end  
                if isnumeric(maxSpace)
                    thisField.maxSpace=maxSpace;

                else
                    error('Maximum smapled radius must be numeric')
                end  
                if isnumeric(FWHM)
                    thisField.FWHM=FWHM;
                else
                    error('FWHM must be numeric')
                end  
                if isnumeric(width)
                    thisField.Width=width;
                else
                    error('Width must be numeric')
                end
                if isnumeric(energy)
                    thisField.Energy=energy;
                else
                    error('Energy must be numeric')
                end
                if isnumeric(dispersionCoefficients)
                    thisField.DispersionCoefficients=dispersionCoefficients;
                else
                    error('DispersionCoefficients must be numeric')
                end
                F_maximalni_necentrovana = 2.9979e8/minWavelength;  %cim vetsi F max tim mensi casovy vektor
                F_maximalni = F_maximalni_necentrovana-thisField.CentralFrequency;
                DF = 2*(F_maximalni);
                thisField.FrequencyVector=linspace(-DF/2,DF/2,dimension(2));
                Dt = (dimension(2)-1)/DF;
                thisField.TimeVector=linspace(-Dt/2,Dt/2,dimension(2));
                dt=Dt/dimension(2);
                f_shifted=thisField.FrequencyVector+thisField.CentralFrequency;
                thisField.WavelengthVector=2.9979e8./f_shifted;
                R = thisField.maxSpace;          % Maximum sampled radius
                N = dimension(1);        % Number of sampling points
                ord = 0;        % Transformation order
                c=load('c_long');
                c=c.c3;
                c = c(ord+1,1:N+1);
                V = c(N+1)/(2*pi*R);    % Maximum frequency
                thisField.SpaceVector = c(1:N)'*R/c(N+1);   % Radius vector                
                thisField.SpaceFrequencyVector = c(1:N)'/(2*pi*R);   % Frequency vector
                r_display = linspace(-max(thisField.SpaceVector),max(thisField.SpaceVector),2*length(thisField.SpaceVector)-1);
                [thisField.MeshSpace,thisField.MeshTime]= meshgrid(r_display,thisField.TimeVector);
                [Jn,Jm] = meshgrid(c(1:N),c(1:N));
                thisField.HenkelC = (2/c(N+1))*besselj(ord,Jn.*Jm/c(N+1))./(abs(besselj(ord+1,Jn)).*abs(besselj(ord+1,Jm)));
                thisField.HenkelInputTransform = (abs(besselj(ord+1,c(1:N)))/R)';   %% prepares input vector for transformation
                thisField.HenkelOuputTransform = thisField.HenkelInputTransform*R/V;  %% prepares output vector for display
                Er = exp(-(thisField.SpaceVector).^2/width^2);       
                Amplitude = sqrt(4*energy*sqrt(1.386*2)/(2.6531*1e-3*FWHM*sqrt(pi)*pi*width^2));
                Et =  exp(-1.386*(thisField.TimeVector/FWHM).^2);
                data = zeros(dimension(1),dimension(2));
                for l = 1:dimension(1)
                     for j = 1:dimension(2)
                       data(l,j) = Amplitude*Er(l)*Et(j);
                     end    
                end
                total_energy = 2.6531*1e-3*1/2*sum((thisField.SpaceVector(2)-thisField.SpaceVector(1))*2*pi*thisField.SpaceVector.*sum(abs(data).^2*(thisField.TimeVector(2)-thisField.TimeVector(1)),2));
                data = sqrt(thisField.Energy)*data/sqrt(total_energy);
                if dispersionCoefficients(1) ~= 0 || dispersionCoefficients(2) ~= 0 || dispersionCoefficients(3) ~= 0
                    A1=thisField.Size(2)*ifftshift(dt*ifft(ifftshift(data,2),thisField.Size(2),2),2);       
                    R=exp(1i*thisField.DispersionCoefficients(1)*((thisField.FrequencyVector*2*pi).^2)/2);
                    R2=exp(1i*thisField.DispersionCoefficients(2)*((thisField.FrequencyVector*2*pi).^3)/6);
                    R3=exp(1i*dispersionCoefficients(3)*((thisField.FrequencyVector*2*pi).^4)/24);
                    for j=1:thisField.Size(1)
                        A1(j,:) =  A1(j,:).*R.*R2;  
                    end  
                    a1 = 1/thisField.Size(2)*fftshift(DF*fft(fftshift(A1,2),thisField.Size(2),2),2);           
                    thisField.InitialSpectrum=A1(1,:);
                    thisField.Data=a1;       
                else
                    thisField.Data=data;
                    thisField.InitialSpectrum=ifftshift(ifft(ifftshift(data(1,:))));
                end    
                thisField.bIntegral=0;
                
            elseif length(fieldnames(InitialParameters))~=8
                error('Incorrect ammount of input arguments')
            end    
        end
        function temporalEnvelope = get.TemporalEnvelope(ThisField)
            temporalEnvelope=abs(ThisField.Data(1,:)).^2;
        end
        function spaceEnvelope = get.SpaceEnvelope(ThisField)
            spaceEnvelope=abs(ThisField.Data(:,round(ThisField.Size(2)/2))).^2;
        end    
        function frequencySpectrum = get.FrequencySpectrum(ThisField)
            frequencySpectrum = ifftshift(ifft(ifftshift(ThisField.Data(1,:))));
        end
        function spaceFrequencySpectrum = get.SpaceFrequencySpectrum(ThisField)
            for i=1:ThisField.Size(2)            
                A1(:,i) = ThisField.Data(:,i)./ThisField.HenkelInputTransform;     %% Prepare vector for transformation
                A1(:,i) = ThisField.HenkelC*A1(:,i);                   %% Obtain the Hankel transform
                A1(:,i) = A1(:,i).*ThisField.HenkelOuputTransform;                 %% Prepare vector for display
            end
            spaceFrequencySpectrum=A1(:,round(ThisField.Size(2)/2));
        end
        function plotTime(ThisField,Element,figureHandle)
            plotHandle=figureHandle;
            plot(ThisField.TimeVector*1e15,1/2*Element.RefractiveIndex*2.6531*1e-3*abs(ThisField.Data(1,:)).^2*1e-9*1e-4);
%             plot(ThisField.TimeVector*1e15,real(ThisField.Data(1,:))/max(real(ThisField.Data(1,:)))*max(1/2*Element.RefractiveIndex*2.65*1e-3*ThisField.TemporalEnvelope*1e-9*1e-4))
            set(plotHandle, 'xlim', [-2*ThisField.FWHM*1e15 2*ThisField.FWHM*1e15]);
%             set(plotHandle, 'xlim', [-3000 3000]);
%             set(plotHandle, 'ylim', [1 3e4]);
           titulek = 'Temporal envelope';
            title(titulek);
            xlabel('Time [fs]')
            ylabel('Intensity [GW/cm^2]')
            M=max(ThisField.TemporalEnvelope);  
            G2 = find(ThisField.TemporalEnvelope>M/10);             %tim zjistime kde priblizne budeme interpolovat
            tq = linspace(ThisField.TimeVector(min(G2)),ThisField.TimeVector(max(G2)),5000);
            attq = interp1(ThisField.TimeVector,ThisField.TemporalEnvelope,tq);
            G3 = find(attq>max(attq)/2);
            fwhm = (max(G3)-min(G3))*(tq(2)-tq(1));
            ThisField.FWHM=fwhm;
            titulek = sprintf('Temporal Envelope FWHM=%.3g fs',fwhm*1e15);
            xlim([-5*fwhm*1e15,5*fwhm*1e15])
            title(titulek);
        end
        function plotTimeSimple(ThisField,RefractiveIndex,figureHandle)
            plotHandle=figureHandle;
            plot(ThisField.TimeVector*1e15,1/2*RefractiveIndex*2.6531*1e-3*ThisField.TemporalEnvelope*1e-9*1e-4);
%             set(plotHandle, 'xlim', [-45*ThisField.FWHM*1e15 45*ThisField.FWHM*1e15]);
            set(plotHandle, 'xlim', [-3e12 3e12]);
            titulek = 'Temporal envelope';
            title(titulek);
            xlabel('Time [fs]')
            ylabel('Intensity [GW/cm^2]')
            M=max(ThisField.TemporalEnvelope);  
            G2 = find(ThisField.TemporalEnvelope>M/10);             %tim zjistime kde priblizne budeme interpolovat
            tq = linspace(ThisField.TimeVector(min(G2)),ThisField.TimeVector(max(G2)),5000);
            attq = interp1(ThisField.TimeVector,ThisField.TemporalEnvelope,tq);
            G3 = find(attq>max(attq)/2);
            fwhm = (max(G3)-min(G3))*(tq(2)-tq(1));
            ThisField.FWHM=fwhm;
            titulek = sprintf('Temporal Envelope FWHM=%.3g fs',fwhm*1e15);
%             xlim([-5*fwhm*1e15,5*fwhm*1e15])
            title(titulek);
        end
        function index = findIndex(ThisField,data,value)
            [minimum,index] = min(abs(data-value));           
        end
        function newData=smoothEdge(ThisField,left,right,data,dimension)
            x11=1:1:left;
            wf1=-(-1-tanh(-5+10*x11/(left)))/2;
            x22=1:1:(right);
            wf2=(1-tanh(-5+10*x22/(right)))/2;
            newData=data;
            newData(:,1:left)=wf1.*newData(:,1:left);
            newData(:,(size(data,dimension)-right+1):end)=wf2.*newData(:,(size(data,dimension)-right+1):end);

        end
        function plotSpace(ThisField,Element,figureHandle)
            plotHandle=figureHandle;
            plot(abs(ThisField.Data(:,ThisField.Size(2)/2+1)).^2)
%             x=linspace(-max(ThisField.SpaceVector)/2,max(ThisField.SpaceVector)/2,200);
%             y=linspace(-max(ThisField.SpaceVector)/2,max(ThisField.SpaceVector)/2,200);
%             [X,Y]= meshgrid(x,y);
%             A=zeros(length(x),length(x));
%             for i=1:length(x)
%                 for j=1:length(x)
%                     r=sqrt(X(i,j)^2+Y(i,j)^2);
%                     index=findIndex(ThisField,ThisField.SpaceVector,r);
%                     A(i,j)=abs(ThisField.Data(index,round(ThisField.Size(2)/2))).^2;
%                 end
%             end    
%             A = squeeze(A(:,:));
%             h=surf(X,Y,A);
%             view(0,90)
%             set(h,'LineStyle','none')
%             titulek = sprintf('Spatial Profile');
%             title(titulek);
%             xlabel('x [m]')
%             ylabel('y [m]')
        end    
        function plotRadius(ThisField,figureHandle,ZAxis)
            rq=linspace(min(ThisField.SpaceVector),max(ThisField.SpaceVector),5e4);
            intensityInSpaceInterval=zeros(1,ThisField.Size(1));
            for j=1:ThisField.Size(1)-1
%                 intensityInSpaceInterval(j)=ThisField.SpaceVector(j+1)*(ThisField.SpaceVector(j+1)-ThisField.SpaceVector(j))*2*pi*sum((ThisField.TimeVector(2)-ThisField.TimeVector(1))*abs(ThisField.Data(j+1,:)).^2);                
                intensityInSpaceInterval(j)=sum((ThisField.TimeVector(2)-ThisField.TimeVector(1))*abs(ThisField.Data(j+1,:)).^2); 
%                 intensityInSpaceInterval(j)=abs(ThisField.Data(j+1,2^13/2)).^2; 
  
            end        
            aq=interp1(ThisField.SpaceVector,1/2*2.6531*1e-3*intensityInSpaceInterval,rq);
            M=max(aq(:));  
            G=find(aq(:)>M/7.389);     
            G2=find(aq(:)>M/2);
            G3=find(aq(:)>M/100); 
            width=rq(max(G));
            width2=rq(max(G2));
            
            % t indep case
            r = ThisField.SpaceVector; %starts at zero, faster to just interp first quadrant then concatenate
            input = intensityInSpaceInterval;
            % make array 
            M = sqrt(r.^2+r'.^2); %dummy matrix of different radii from origin
            M = reshape(M, [1,length(r)^2]); %make into 1D list
            output = interp1(r, input, M,'spline',0); %interpolate onto 1D list
            output = reshape(output, [length(r),length(r)]); %make back into array
            %put into all four quadrants
            output = [rot90(output); output];
            output = [rot90(rot90(output)) output];
            r_full = linspace(-max(r),max(r),2*length(r));
            width3=2*sqrt(sum(sum((ThisField.SpaceVector(2)-ThisField.SpaceVector(1))^2*r_full.^2.*output))/sum(sum((ThisField.SpaceVector(2)-ThisField.SpaceVector(1))^2.*output)));
            ThisField.WidthVector(nnz(ThisField.WidthVector)+1) = width;
            ThisField.WidthVector2(nnz(ThisField.WidthVector)+1) = width2;
            ThisField.WidthVector3(nnz(ThisField.WidthVector)+1) = width3;
            plotHandle=figureHandle;
            plot(ZAxis(1:nnz(ZAxis))*1000,1000*ThisField.WidthVector(1:nnz(ZAxis)),ZAxis(1:nnz(ZAxis))*1000,1000*ThisField.WidthVector2(1:nnz(ZAxis)),ZAxis(1:nnz(ZAxis))*1000,1000*ThisField.WidthVector3(1:nnz(ZAxis)))

            set(plotHandle, 'xlim', [ZAxis(1)*1000 (2*max(ZAxis)+0.001)*1000]);
%             set(plotHandle, 'ylim', [0 0.1]);
            titulek = sprintf('Beam radius integrated=%.3g mm Propagated Distance=%.3g mm',width*1000,max(ZAxis)*1000 );
            title(titulek);
            xlabel('Propagated distance [mm]')
            ylabel('Beam radius 1/e^2 [mm]')
            legend('1/e^2','FWHM','D4sigma')
            drawnow;
        end   
        function plotSpaceTime(ThisField,Element,figureHandle)
            plotHandle=figureHandle;
            A_disp = zeros(2*ThisField.Size(1)-1,ThisField.Size(2)); 
            A_disp(ThisField.Size(1):end,:)= ThisField.Data;
            for i=1:ThisField.Size(1)-1
                A_disp(i,:) = ThisField.Data(end-i,:);
            end
            intensityXY = squeeze(A_disp(:,:));
            A=abs(intensityXY').^2/max(max(abs(intensityXY').^2));
%             A=A+1e-4;
%             h=mesh(ThisField.MeshSpace,1*1e12*ThisField.MeshTime,A,log10(A));
            h=mesh(ThisField.MeshSpace,1*1e12*ThisField.MeshTime,A);
            view(0,90)
%             set(gca, 'ZScale', 'log');
            titulek = sprintf('Pulse Energy %.3g mJ; B integral %.2g',ThisField.Energy*1000, ThisField.bIntegral);
            title(titulek);
            xlabel('x [m]')
            ylabel('t [ps]')
%             set(plotHandle, 'xlim', [-0.025 0.025]);
%             set(plotHandle, 'ylim', [-3 3]);
%             cb2 = colorbar;
%             axpos = h.Parent.Position;
%             axes;
%             axis off
%             cb = colorbar('Position',cb2.Position);
%             caxis(log10([min(A(A>0)) max(A(:))]))
%             cb.TickLabels = sprintf('10^{%1.1f}\n',cb.Ticks);
%             delete(cb2)
%             set(h.Parent,'Position',axpos)
        end
        function plotSpectrum(ThisField,figureHandle)
            plotHandle=figureHandle;     
            plot(ThisField.WavelengthVector*1e9,ThisField.FrequencySpectrum.*conj(ThisField.FrequencySpectrum).*(2*pi*(ThisField.FrequencyVector+ThisField.CentralFrequency)).^2/max(ThisField.FrequencySpectrum.*conj(ThisField.FrequencySpectrum).*(2*pi*(ThisField.FrequencyVector+ThisField.CentralFrequency)).^2),ThisField.WavelengthVector*1e9,ThisField.InitialSpectrum.*conj(ThisField.InitialSpectrum)/max(ThisField.InitialSpectrum.*conj(ThisField.InitialSpectrum)));
%             plot(ThisField.WavelengthVector*1e9,ThisField.AverageSpectrumPropagate(nnz(ThisField.ZAxis),:)/max(ThisField.AverageSpectrumPropagate(nnz(ThisField.ZAxis),:)),ThisField.WavelengthVector*1e9,ThisField.InitialSpectrum.*conj(ThisField.InitialSpectrum)/max(ThisField.InitialSpectrum.*conj(ThisField.InitialSpectrum)));
%         plot(ThisField.WavelengthVector*1e9,ThisField.FrequencySpectrum.*conj(ThisField.FrequencySpectrum).*(2*pi*(ThisField.FrequencyVector+ThisField.CentralFrequency)).^2/max(ThisField.FrequencySpectrum.*conj(ThisField.FrequencySpectrum).*(2*pi*(ThisField.FrequencyVector+ThisField.CentralFrequency)).^2));

            
%             set(plotHandle, 'xlim', [ThisField.MinWavelength*1e9 (2*ThisField.CentralWavelength-ThisField.MinWavelength)*1e9]);
            set(plotHandle, 'xlim', [700 1000]);

%             set(plotHandle, 'ylim', [1e-7 1.1]);

            titulek3 = 'Power spectrum at r = 0 (initial in red, output in blue)' ;
            title(titulek3);
            xlabel('Wavelength [nm]')
            ylabel('A.U.')

        end                  
        function plotSpaceSpectrum(ThisField,figureHandle)
            plotHandle=figureHandle; 
            A1 = ThisField.Data;
%             A1 = ifftshift((ThisField.TimeVector(2)-ThisField.TimeVector(1))*ifft(((ifftshift(ThisField.Data(:,2^10/2))))));
            InputMatrix=repmat(ThisField.HenkelInputTransform,1,ThisField.Size(2));
            OutputMatrix=repmat(ThisField.HenkelOuputTransform,1,ThisField.Size(2));
            A1 = 1./InputMatrix.*A1;     %% Prepare vector for transformation
            A1 = ThisField.HenkelC*A1;                   %% Obtain the Hankel transform
            A1 = OutputMatrix.*A1;   
            plot(abs(A1(:,round(ThisField.Size(2)/2))).^2);
%             set(plotHandle, 'xlim', [0 1000]);
            titulek3 = 'Space spectrum' ;
            xlabel('k-vectors')
            title(titulek3);
        end    
        function plotSpaceSpectrumAngle(ThisField,figureHandle)
            plotHandle=figureHandle;            
            for i=1:ThisField.Size(1)
                Wavefront(i,:)=fftshift(fft(fftshift(ThisField.Data(i,:))));
            end             
            plot(1e3*ThisField.SpaceVector,unwrap(angle(Wavefront(:,round(ThisField.Size(2)/2)))));
%             set(plotHandle, 'xlim', [0 0.02]);
            titulek3 = 'Wavefront' ;
            title(titulek3);
            xlabel('x [mm]')
            ylabel('Angle [rad]')
        end  
        function changeSpaceSampling(ThisField,newMaxSpace)
            ThisField.maxSpace=newMaxSpace.maxR;
            oldR=ThisField.SpaceVector;
            R = ThisField.maxSpace;          % Maximum sampled radius
            N = ThisField.Size(1);        % Number of sampling points
            ord = 0;        % Transformation order
            c=load('c_long');
            c=c.c3;
            c = c(ord+1,1:N+1);
            V = c(N+1)/(2*pi*R);    % Maximum frequency
            ThisField.SpaceVector = c(1:N)'*R/c(N+1);   % Radius vector                
            ThisField.SpaceFrequencyVector = c(1:N)'/(2*pi*R);   % Frequency vector
            r_display = linspace(-max(ThisField.SpaceVector),max(ThisField.SpaceVector),2*length(ThisField.SpaceVector)-1);
            [ThisField.MeshSpace,ThisField.MeshTime]= meshgrid(r_display,ThisField.TimeVector);
            [Jn,Jm] = meshgrid(c(1:N),c(1:N));
            ThisField.HenkelC = (2/c(N+1))*besselj(ord,Jn.*Jm/c(N+1))./(abs(besselj(ord+1,Jn)).*abs(besselj(ord+1,Jm)));
            % C is the transformation matrix
            ThisField.HenkelInputTransform = (abs(besselj(ord+1,c(1:N)))/R)';   %% prepares input vector for transformation
            ThisField.HenkelOuputTransform = ThisField.HenkelInputTransform*R/V;  %% prepares output vector for display
            interpolationPhase=interp1(oldR,unwrap(angle(ThisField.Data)),ThisField.SpaceVector,'pchip');
            interpolationAmpl=interp1(oldR,abs(ThisField.Data),ThisField.SpaceVector,'pchip');
            interpolation=interpolationAmpl.*exp(1i*interpolationPhase);
%             ThisField.Data=fillmissing(interpolation,'spline');
            ThisField.Data=interpolation;
        end    
        function changeTimeInterval(ThisField,N)
            oldTime=ThisField.TimeVector;
%             minWavelength=min(ThisField.WavelengthVector);
            minWavelength=350e-9;
            F_maximalni_necentrovana = 2.9979e8/minWavelength;  %cim vetsi F max tim mensi casovy vektor
            F_maximalni = F_maximalni_necentrovana-ThisField.CentralFrequency;
            DF = 2*(F_maximalni);
            ThisField.FrequencyVector=linspace(-DF/2,DF/2,N);
            Dt = (N-1)/DF;
            ThisField.TimeVector=linspace(-Dt/2,Dt/2,N);
            dt=Dt/N;
            f_shifted=ThisField.FrequencyVector+ThisField.CentralFrequency;
            ThisField.WavelengthVector=2.9979e8./f_shifted;
            [~,index]=find(ThisField.WavelengthVector>5e-6);
            ThisField.WavelengthVector(1:max(index))=5e-6;
            ThisField.IndexSpectralCutoff(1)=max(index);
            [~,index]=find(ThisField.WavelengthVector<0.490e-6);
            ThisField.IndexSpectralCutoff(2)=min(index);
            ThisField.Size(2)=N;
            for ii=1:ThisField.Size(1)
            interpolation(ii,:)=interp1(oldTime,ThisField.Data(ii,:),ThisField.TimeVector,'pchip');
            end
            interpolation(:,ThisField.TimeVector>max(oldTime))=0*interpolation(:,ThisField.TimeVector>max(oldTime));
            interpolation(:,ThisField.TimeVector<min(oldTime))=0*interpolation(:,ThisField.TimeVector<min(oldTime));
           
            ThisField.Data=fillmissing(interpolation,'pchip');
        end    
        function simulateFocalSpot(ThisField,f)
            focus=f; %%focus of the lens, focusScan je jakoby druha cocka predu focusem ktera focus zkracuje/prodluzuje
            A1 = zeros(ThisField.Size);
%             a1 = zeros(ThisField.Size);
            dt=ThisField.TimeVector(2)-ThisField.TimeVector(1);
            DF=max(ThisField.FrequencyVector)-min(ThisField.FrequencyVector);
            A1=ThisField.Size(2)*ifftshift(dt*ifft(ifftshift(ThisField.Data,2),ThisField.Size(2),2),2);       
            InputMatrix=repmat(ThisField.HenkelInputTransform,1,ThisField.Size(2));
            OutputMatrix=repmat(ThisField.HenkelOuputTransform,1,ThisField.Size(2));
            %             for i=1:ThisField.Size(2)            
            A1 = 1./InputMatrix.*A1;     %% Prepare vector for transformation
            A1 = ThisField.HenkelC*A1;                   %% Obtain the Hankel transform
            A1 = OutputMatrix.*A1;                 %% Prepare vector for display
            for k =1:ThisField.Size(2)
                OutputSpaceVector{k}=ThisField.SpaceFrequencyVector*(ThisField.WavelengthVector(k)*focus);
            end        
            for m=1:ThisField.Size(2)
                for k=1:ThisField.Size(1)
                    A1(k,m) = A1(k,m).*exp(1i*(pi/(ThisField.WavelengthVector(m)*focus))*(OutputSpaceVector{m}(k)).^2*focus^2*ThisField.WavelengthVector(m)^2)./(1i*ThisField.WavelengthVector(m)*focus);  %fourierka pro t v kazdem bode u,v
                end
            end
        
            PeakEnvelopeFocus = ifftshift(ifft(((ifftshift(A1(1,:))))));  %fourierka pro t v kazdem bode u,v
            WholeFocus=1/ThisField.Size(2)*fftshift(DF*fft(fftshift(A1,2),ThisField.Size(2),2),2);
            ThisField.Data=WholeFocus;
            %zmena space vektoru
            R = max(OutputSpaceVector{round(ThisField.Size(2)/2)});         % Maximum sampled radius
            N = ThisField.Size(1);        % Number of sampling points
            ord = 0;        % Transformation order
            c=load('c_long');
            c=c.c3;
            c = c(ord+1,1:N+1);
            V = c(N+1)/(2*pi*R);    % Maximum frequency
            ThisField.SpaceVector = c(1:N)'*R/c(N+1);   % Radius vector                
            ThisField.SpaceFrequencyVector = c(1:N)'/(2*pi*R);   % Frequency vector
            r_display = linspace(-max(ThisField.SpaceVector),max(ThisField.SpaceVector),2*length(ThisField.SpaceVector)-1);
            [ThisField.MeshSpace,ThisField.MeshTime]= meshgrid(r_display,ThisField.TimeVector);
            [Jn,Jm] = meshgrid(c(1:N),c(1:N));
            ThisField.HenkelC = (2/c(N+1))*besselj(ord,Jn.*Jm/c(N+1))./(abs(besselj(ord+1,Jn)).*abs(besselj(ord+1,Jm)));
            % C is the transformation matrix
            ThisField.HenkelInputTransform = (abs(besselj(ord+1,c(1:N)))/R)';   %% prepares input vector for transformation
            ThisField.HenkelOuputTransform = ThisField.HenkelInputTransform*R/V; 
        end    
        function simulateFocalSpotAfterLens(ThisField,f,distanceFromLens)
            focus=f; %%focus of the lens, focusScan je jakoby druha cocka predu focusem ktera focus zkracuje/prodluzuje
            A1 = zeros(ThisField.Size);
%             a1 = zeros(ThisField.Size);
            dt=ThisField.TimeVector(2)-ThisField.TimeVector(1);
            DF=max(ThisField.FrequencyVector)-min(ThisField.FrequencyVector);
            A1=ThisField.Size(2)*ifftshift(dt*ifft(ifftshift(ThisField.Data,2),ThisField.Size(2),2),2);       
            InputMatrix=repmat(ThisField.HenkelInputTransform,1,ThisField.Size(2));
            OutputMatrix=repmat(ThisField.HenkelOuputTransform,1,ThisField.Size(2));
            %             for i=1:ThisField.Size(2)            
            A1 = 1./InputMatrix.*A1;     %% Prepare vector for transformation
            A1 = ThisField.HenkelC*A1;                   %% Obtain the Hankel transform
            A1 = OutputMatrix.*A1;                 %% Prepare vector for display
            for k =1:ThisField.Size(2)
                OutputSpaceVector{k}=ThisField.SpaceFrequencyVector*(ThisField.WavelengthVector(k)*focus);
            end        
            for m=1:ThisField.Size(2)
                for k=1:ThisField.Size(1)
                    A1(k,m) = A1(k,m).*exp(1i*(1-distanceFromLens/focus)*(pi/(ThisField.WavelengthVector(m)*focus))*(OutputSpaceVector{m}(k)).^2*focus^2*ThisField.WavelengthVector(m)^2)./(1i*ThisField.WavelengthVector(m)*focus);  %fourierka pro t v kazdem bode u,v
                end
            end
        
            PeakEnvelopeFocus = ifftshift(ifft(((ifftshift(A1(1,:))))));  %fourierka pro t v kazdem bode u,v
            WholeFocus=1/ThisField.Size(2)*fftshift(DF*fft(fftshift(A1,2),ThisField.Size(2),2),2);
            ThisField.Data=WholeFocus;
            %zmena space vektoru
            R = max(OutputSpaceVector{round(ThisField.Size(2)/2)});         % Maximum sampled radius
            N = ThisField.Size(1);        % Number of sampling points
            ord = 0;        % Transformation order
            c=load('c_long');
            c=c.c3;
            c = c(ord+1,1:N+1);
            V = c(N+1)/(2*pi*R);    % Maximum frequency
            ThisField.SpaceVector = c(1:N)'*R/c(N+1);   % Radius vector                
            ThisField.SpaceFrequencyVector = c(1:N)'/(2*pi*R);   % Frequency vector
            r_display = linspace(-max(ThisField.SpaceVector),max(ThisField.SpaceVector),2*length(ThisField.SpaceVector)-1);
            [ThisField.MeshSpace,ThisField.MeshTime]= meshgrid(r_display,ThisField.TimeVector);
            [Jn,Jm] = meshgrid(c(1:N),c(1:N));
            ThisField.HenkelC = (2/c(N+1))*besselj(ord,Jn.*Jm/c(N+1))./(abs(besselj(ord+1,Jn)).*abs(besselj(ord+1,Jm)));
            % C is the transformation matrix
            ThisField.HenkelInputTransform = (abs(besselj(ord+1,c(1:N)))/R)';   %% prepares input vector for transformation
            ThisField.HenkelOuputTransform = ThisField.HenkelInputTransform*R/V; 
        end    
        function linearStep(ThisField,Element)    
            
            c_0=2.9979e8;
            dt=ThisField.TimeVector(2)-ThisField.TimeVector(1);
            DF=max(ThisField.FrequencyVector)-min(ThisField.FrequencyVector);
            A1=ThisField.Size(2)*ifftshift(dt*ifft(ifftshift(ThisField.Data,2),ThisField.Size(2),2),2);       
            InputMatrix=repmat(ThisField.HenkelInputTransform,1,ThisField.Size(2));
            OutputMatrix=repmat(ThisField.HenkelOuputTransform,1,ThisField.Size(2));
            %             for i=1:ThisField.Size(2)            
            A1 = 1./InputMatrix.*A1;     %% Prepare vector for transformation
            A1 = ThisField.HenkelC*A1;                   %% Obtain the Hankel transform
            A1 = OutputMatrix.*A1;                 %% Prepare vector for display
%             end

%             zkouska=A1;
            f_shifted=ThisField.FrequencyVector+ThisField.CentralFrequency;    
            coeffs=getDispersionCoefficients(Element);
            k1=Element.RefractiveIndex*ThisField.WaveNumber+coeffs(1)/Element.Length*(2*pi*ThisField.FrequencyVector);
            k=(2*pi*f_shifted.*Element.RefractiveIndexVector./c_0);
%             k(1:2200)=k(2200);
%             k(6000:end)=k(6000);
%             k1(1:2200)=k1(2200);
%             k1(6000:end)=k1(6000);
%             k1(1:ThisField.IndexSpectralCutoff(1))=k1(ThisField.IndexSpectralCutoff(1));
%             k1(ThisField.IndexSpectralCutoff(2):end)=k1(ThisField.IndexSpectralCutoff(2));
%             k(1:ThisField.IndexSpectralCutoff(1))=k(ThisField.IndexSpectralCutoff(1));
%             k(ThisField.IndexSpectralCutoff(2):end)=k(ThisField.IndexSpectralCutoff(2));
% % k=(2*pi*Element.RefractiveIndexVector./ThisField.WavelengthVector);
            for m=1:ThisField.Size(2)
%                   if m==2048
%                       m=m;
%                   end    
%                       
                  odm=k(m)^2-4*pi^2*ThisField.SpaceFrequencyVector.^2;
                  odm(odm<0)=0;
% %               
                  H = exp(1i*Element.Dz*sqrt(odm));
                  H2 = exp(-1i*Element.Dz*k1(m));
                  A1(:,m) = A1(:,m).*H.*H2;    
%                   
                  
%                   H = exp(1i*(k(m)^2-k1(m)^2)*Element.Dz/(2*k1(m)));
%                   H2 = exp(-1i*4*pi^2*Element.Dz*ThisField.SpaceFrequencyVector.^2./(2*k1(m)));
%                   A1(:,m) = A1(:,m).*H.*H2;   


            end

%             for i=1:ThisField.Size(2) 

%             a1 = zeros(ThisField.Size);
            a1 = ThisField.HenkelC*(1./OutputMatrix.*A1);                    %% Obtain the Hankel transform
            a1 = InputMatrix.*a1; %% Prepare vector for display
%             end
%             x1=linspace(50,0,ThisField.IndexSpectralCutoff(1));
%             x2=linspace(0,50,ThisField.Size(2)-ThisField.IndexSpectralCutoff(2)+1);
%             x11=1:1:ThisField.IndexSpectralCutoff(1);
%             wf1=-(-1-tanh(-5+10*x11/(ThisField.IndexSpectralCutoff(1))))/2;
%             x22=1:1:(ThisField.Size(2)-ThisField.IndexSpectralCutoff(2)+1);
%             wf2=(1-tanh(-5+10*x22/(ThisField.Size(2)-ThisField.IndexSpectralCutoff(2)+1)))/2;
%             a1(:,1:ThisField.IndexSpectralCutoff(1))=wf1.*a1(:,1:ThisField.IndexSpectralCutoff(1));
% %             a1(:,ThisField.IndexSpectralCutoff(2):end)=exp(-x2).*a1(:,ThisField.IndexSpectralCutoff(2):end);
%             a1(:,ThisField.IndexSpectralCutoff(2):end)=wf2.*a1(:,ThisField.IndexSpectralCutoff(2):end);
% %             r1=linspace(0,20,11);
%             x=1:1:ThisField.Size(1)/8;
% %             w=0.5*(1-cos(2*pi*x/ThisField.Size(1))); Hamming window
%             w=(1-tanh(-5+10*x/(ThisField.Size(1)/8)))/2;
%             a1(end-ThisField.Size(1)/8+1:end,:)=transpose(w).*a1(end-ThisField.Size(1)/8+1:end,:);
%             
%       
            
            
            ThisField.Data = 1/ThisField.Size(2)*fftshift(DF*fft(fftshift(a1,2),ThisField.Size(2),2),2);
%             x11=1:1:40;
%             wf1=-(-1-tanh(-5+10*x11/(40)))/2;
%             x22=1:1:(40);
%             wf2=(1-tanh(-5+10*x22/(40)))/2;
%             ThisField.Data(:,1:40)=wf1.*ThisField.Data(:,1:40);
%             ThisField.Data(:,(ThisField.Size(2)-40+1):end)=wf2.*ThisField.Data(:,(ThisField.Size(2)-40+1):end);
% 
%             
%             
%             
            
            total_energy_time = Element.RefractiveIndex*2.6531*1e-3*1/2*sum((ThisField.SpaceVector(2)-ThisField.SpaceVector(1))*2*pi*ThisField.SpaceVector.*sum(abs(ThisField.Data).^2*(ThisField.TimeVector(2)-ThisField.TimeVector(1)),2));
            %calculate energy and b intergral
            ThisField.Energy=total_energy_time;
            ThisField.bIntegral=ThisField.bIntegral+2*pi/ThisField.CentralWavelength*Element.NonlinearIndex*Element.RefractiveIndex*2.6531*1e-3*1/2*max(max(abs(ThisField.Data).^2))*Element.Dz;
            %calculate spectrum and average spectrum
            
            
            average_spectrum = zeros(1,ThisField.Size(2));
            r=ThisField.SpaceVector;            
            average_spectrum = average_spectrum+Element.RefractiveIndex*(r(1))/r(end)*2*pi*r(1)*abs(a1(1,:)).^2;
            for i =2:ThisField.Size(1)
            average_spectrum = average_spectrum+Element.RefractiveIndex*(r(i)-r(i-1))/r(end)*2*pi*r(i)*abs(a1(i,:)).^2;
            end
%             for i =2:10
%             average_spectrum2 = average_spectrum+(r(i)-r(i-1))/r(end)*2*pi*r(i)*abs(a1(i,:)).^2;
%             end
            average_spectrum = average_spectrum.*(2*pi*(ThisField.FrequencyVector+ThisField.CentralFrequency)).^2;
%             average_spectrum2 = average_spectrum2/10.*(2*pi*(ThisField.FrequencyVector+ThisField.CentralFrequency)).^2;
%             ThisField.SpectrumPropagate(nnz(ThisField.ZAxis),:)=a1(1,:);
            ThisField.AverageSpectrumPropagate(nnz(ThisField.ZAxis),:)=average_spectrum;          
            ThisField.TimePropagate(nnz(ThisField.ZAxis),:)=ThisField.Data(1,:);
        end                        
        function nonlinearStep(ThisField,Element)          
            dt=ThisField.TimeVector(2)-ThisField.TimeVector(1);
            DF=max(ThisField.FrequencyVector)-min(ThisField.FrequencyVector);
            Dz=Element.Dz;
            RungeKuttaSteps=Element.RungeKuttaSteps;
            n1=Element.RefractiveIndex;
            n2=Element.NonlinearIndex;
            L_0=ThisField.CentralWavelength;
            k_0=2*pi/ThisField.CentralWavelength;
            w_0=2*pi*ThisField.CentralFrequency;
            c_0=2.9979e8;          
            f=ThisField.FrequencyVector;  
            f_shifted=c_0./ThisField.WavelengthVector;
            dimension=ThisField.Size;
            kv=2*pi./ThisField.WavelengthVector;
            z=linspace(0,Dz,RungeKuttaSteps);
            h=z(2)-z(1);
        A_z=dimension(2)*ifftshift(dt*ifft(ifftshift(ThisField.Data,2),ThisField.Size(2),2),2);
        
        coeffs=getDispersionCoefficients(Element);
        k1=Element.RefractiveIndex*ThisField.WaveNumber+coeffs(1)/Element.Length*(2*pi*ThisField.FrequencyVector);
        maximumTime=zeros(1,ThisField.Size(1));
        for ppp=1:ThisField.Size(1) 
            maximumTime(ppp)=max((n1*2.6531*1e-3/2)*abs(ThisField.Data(ppp,:)).^2);
        end    
        maximumTime=maximumTime*1e-9*1e-4;
        whereToSolve=find(maximumTime>0.01);
        for k = 2:length(z)+1  
            F1 = 1/dimension(2)*fftshift(DF*fft(fftshift(A_z(whereToSolve,:),2),ThisField.Size(2),2),2);
            F2 = (n1*2.6531*1e-3/2)*abs(F1).^2;
            F_final = dimension(2)*ifftshift(dt*ifft(ifftshift(F1.*F2,2),ThisField.Size(2),2),2);
            Ka = 1i.*n1.*kv.^2./k1.*n2.*F_final;
            K1=Ka;
            F1 = 1/dimension(2)*fftshift(DF*fft(fftshift(A_z(whereToSolve,:)+h/2*K1,2),ThisField.Size(2),2),2);
            F2 = (n1*2.6531*1e-3/2)*abs(F1).^2;
            F_final = dimension(2)*ifftshift(dt*ifft(ifftshift(F1.*F2,2),ThisField.Size(2),2),2);
            Ka = 1i.*n1.*kv.^2./k1.*n2.*F_final;
            K2=Ka;
            F1 = 1/dimension(2)*fftshift(DF*fft(fftshift(A_z(whereToSolve,:)+h/2*K2,2),ThisField.Size(2),2),2);
            F2 = (n1*2.6531*1e-3/2)*abs(F1).^2;
            F_final = dimension(2)*ifftshift(dt*ifft(ifftshift(F1.*F2,2),ThisField.Size(2),2),2);
            Ka = 1i.*n1.*kv.^2./k1.*n2.*F_final;
            K3=Ka;
            F1 = 1/dimension(2)*fftshift(DF*fft(fftshift(A_z(whereToSolve,:)+h*K3,2),ThisField.Size(2),2),2);
            F2 = (n1*2.6531*1e-3/2)*abs(F1).^2;
            F_final = dimension(2)*ifftshift(dt*ifft(ifftshift(F1.*F2,2),ThisField.Size(2),2),2);
            Ka = 1i.*n1.*kv.^2./k1.*n2.*F_final;
            K4=Ka;
            A_z(whereToSolve,:)= A_z(whereToSolve,:)+h/6*(K1+2*K2+2*K3+K4);
        end

        ThisField.Data = 1/dimension(2)*fftshift(DF*fft(fftshift(A_z,2),ThisField.Size(2),2),2);
        [row, col] = find(isnan((abs(ThisField.Data))));
        ZbavNAN=[row,col];
        if isempty(row)
        else    
            Sizeru=size(ZbavNAN);
            for ii=1:Sizeru(1)
            ThisField.Data(ZbavNAN(ii,1),ZbavNAN(ii,2))=0;
            end
        end    
        
        end
        function propagate(ThisField,OpticalPath,Figure)
            for i=1:length(OpticalPath)
            Element=OpticalPath{i}; % now only for materials,lens,disperion , for mirror needs to be error checked
             
            if isa(Element,'Lens')
                    %the smaller the focus, the bigger must be number of spatial points
                dt=ThisField.TimeVector(2)-ThisField.TimeVector(1);
                DF=max(ThisField.FrequencyVector)-min(ThisField.FrequencyVector);
                A_lens=ThisField.Size(2)*ifftshift(dt*ifft(ifftshift(ThisField.Data,2),ThisField.Size(2),2),2);
                for b = 1:ThisField.Size(1)
%                     for p = 1:ThisField.Size(2)
%                         A_lens(b,p) = A_lens(b,p)*exp((-1i*pi/(Element.Focus*ThisField.WavelengthVector(p)))*ThisField.SpaceVector(b)^2); 
%                     end
                    A_lens(b,:) = A_lens(b,:).*exp((-1i*pi./(Element.Focus*ThisField.WavelengthVector(:)')).*ThisField.SpaceVector(b)^2); 
                end
                ThisField.Data = 1/ThisField.Size(2)*fftshift(DF*fft(fftshift(A_lens,2),ThisField.Size(2),2),2);
                %make sure the energy is conserved
                total_energy = Element.RefractiveIndex*2.6531*1e-3*1/2*sum((ThisField.SpaceVector(2)-ThisField.SpaceVector(1))*2*pi*ThisField.SpaceVector.*sum(abs(ThisField.Data).^2*(ThisField.TimeVector(2)-ThisField.TimeVector(1)),2));
                ThisField.Data = sqrt(ThisField.Energy)*ThisField.Data/sqrt(total_energy);
                %make sure the energy is conserved  
                for j=1:Element.Length/Element.Dz
                    
                    if nnz(ThisField.ZAxis)==0
                        ThisField.ZAxis(nnz(ThisField.ZAxis)+1)=Element.Dz;
                    else    
                        ThisField.ZAxis(nnz(ThisField.ZAxis)+1) = ThisField.ZAxis(nnz(ThisField.ZAxis))+Element.Dz;         
                    end
                    if Element.NonlinearIndex==0
                    else    
                    nonlinearStep(ThisField,Element) 
                    end
                    linearStep(ThisField,Element)
                    figure1=Figure;
                    plothandle4=subplot(3,2,4);
                    plotTime(ThisField,Element,plothandle4)
                    plothandle1=subplot(3,2,1);
                    plotSpaceTime(ThisField,Element,plothandle1)
                    plothandle3=subplot(3,2,3);
                    plotSpectrum(ThisField,plothandle3)
                    plothandle2=subplot(3,2,2);
                    plotRadius(ThisField,plothandle2,ThisField.ZAxis)
                    plothandle5=subplot(3,2,5);
                    plotSpaceSpectrum(ThisField,plothandle5)
                    plothandle6=subplot(3,2,6);
                    plotSpaceSpectrumAngle(ThisField,plothandle6)
             end 
                
             
             elseif isa(Element,'Dispersion')
                    for k =1:ThisField.Size(1)
                        A1(k,:) = ifftshift((ThisField.TimeVector(2)-ThisField.TimeVector(1))*ifft(((ifftshift(ThisField.Data(k,:))))));  %fourierka pro t v kazdem bode u,v
                    end
                    R=exp(1i*Element.GDD*((2*pi*ThisField.FrequencyVector).^2)/2);
                    R2=exp(1i*Element.TOD*((2*pi*ThisField.FrequencyVector).^3)/6);
                    R3=exp(1i*Element.FOD*((2*pi*ThisField.FrequencyVector).^4)/24);
                    for j=1:ThisField.Size(1)
                        A1(j,:) =  A1(j,:).*R.*R2.*R3;  
                    end  
                    DF=max(ThisField.FrequencyVector)-min(ThisField.FrequencyVector);
                    for k =1:ThisField.Size(1)
                        a1(k,:) = fftshift(DF*fft(fftshift(A1(k,:))));  
                    end
                    ThisField.Data=a1; 
            
            elseif isstruct(Element)
                if strcmp(Element.type,'aperture')
                    ThisField.Data(ThisField.SpaceVector>Element.radius,:)=0;
                    ThisField.Energy = 1*2.6531*1e-3*1/2*sum((ThisField.SpaceVector(2)-ThisField.SpaceVector(1))*2*pi*ThisField.SpaceVector.*sum(abs(ThisField.Data).^2*(ThisField.TimeVector(2)-ThisField.TimeVector(1)),2));
                elseif strcmp(Element.type,'changeSpace')
                    changeSpaceSampling(ThisField,Element);
                elseif strcmp(Element.type,'simulateFocus')
                    simulateFocalSpot(ThisField,Element.f);
                elseif strcmp(Element.type,'simulateFocusAfterLens')
                    simulateFocalSpot(ThisField,Element.f);    
                elseif strcmp(Element.type,'changeTime')
                    changeTimeInterval(ThisField,Element.NewN)
                end    
            else    
                            %make sure the energy is conserved
            total_energy = Element.RefractiveIndex*2.6531*1e-3*1/2*sum((ThisField.SpaceVector(2)-ThisField.SpaceVector(1))*2*pi*ThisField.SpaceVector.*sum(abs(ThisField.Data).^2*(ThisField.TimeVector(2)-ThisField.TimeVector(1)),2));
            ThisField.Data = sqrt(ThisField.Energy)*ThisField.Data/sqrt(total_energy);
            %make sure the energy is conserved  
                for j=1:Element.Length/Element.Dz
%                     tic
%                      ThisField.bIntegral
                    if nnz(ThisField.ZAxis)==0
                        ThisField.ZAxis(nnz(ThisField.ZAxis)+1)=Element.Dz;
                    else    
                        ThisField.ZAxis(nnz(ThisField.ZAxis)+1) = ThisField.ZAxis(nnz(ThisField.ZAxis))+Element.Dz;        
                    end
                    if Element.NonlinearIndex==0
                    else    
                    nonlinearStep(ThisField,Element) 
                    end
                    linearStep(ThisField,Element)
                    j
%                     toc
                    figure1=Figure;
                    plothandle4=subplot(3,2,4);
                    plotTime(ThisField,Element,plothandle4)
                    plothandle1=subplot(3,2,1);
                    plotSpaceTime(ThisField,Element,plothandle1)
                    plothandle3=subplot(3,2,3);
                    plotSpectrum(ThisField,plothandle3)
                    plothandle2=subplot(3,2,2);
                    plotRadius(ThisField,plothandle2,ThisField.ZAxis)
                    plothandle5=subplot(3,2,5);
                    plotSpaceSpectrum(ThisField,plothandle5)
                    plothandle6=subplot(3,2,6);
                    plotSpaceSpectrumAngle(ThisField,plothandle6)
                    drawnow;
                end
            end
            end
        end
   end
end    