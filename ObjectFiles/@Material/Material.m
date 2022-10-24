classdef Material < Path
    properties
        Length
        NonlinearIndex
        RefractiveIndex
        RefractiveIndexVector
        WavelengthVector
        Dz
        CentralWavelength
        RungeKuttaSteps
    end


    methods
        function thisMaterial=Material(name,length2,numberOfPoints,rungeKuttaSteps,nonlinearIndex,wavelengthVector,centralWavelength) 
            thisMaterial@Path(name)
            if nargin==7
                if isnumeric(length2)
                    thisMaterial.Length=length2;
                else
                    error('Length must be numeric')
                end
                if isnumeric(rungeKuttaSteps)
                    thisMaterial.RungeKuttaSteps=rungeKuttaSteps;
                else
                    error('Number of steps for Runge Kutta must be numeric')
                end   
                if isnumeric(nonlinearIndex)
                    thisMaterial.NonlinearIndex=nonlinearIndex;
                else
                    error('Nonlinear index must be numeric')
                end  
                if isnumeric(wavelengthVector)
                    thisMaterial.WavelengthVector=wavelengthVector;
                else
                    error('Wavelength vector must be numeric')
                end
                if isnumeric(centralWavelength)
                    thisMaterial.CentralWavelength=centralWavelength;
                else
                    error('Central wavelength must be numeric')
                end 
                if isnumeric(numberOfPoints)
                    thisMaterial.Dz=thisMaterial.Length/numberOfPoints;
                else
                    error('The number of sampling points must be numeric')
                end
                fid =fopen('materials.txt');
                l = sym('l');
                counter=0;
                tline = 'o';
                while ischar(tline)
                    counter = counter + 1;
                    tline = fgets(fid);
                    if tline ~= -1
                        B{counter} = strsplit(tline, ' ', 'CollapseDelimiters',true);
                    else
                    end
                end
                m=0;
                for i=1:length(B)
                    if strcmp(name, B{i}{1})
                        m = i;
                    else
                    end
                end   
                sellmeier = B{m}{2};
                A(l) = evalin(symengine, sellmeier);
                fclose(fid);
                thisMaterial.RefractiveIndexVector = double(A(wavelengthVector));
                thisMaterial.RefractiveIndex = double(A(centralWavelength));
            end
        end
        function  dispersionCoeff = getDispersionCoefficients(Element)
            fid =fopen('materials.txt');
            l = sym('l');
            counter=0;
            tline = 'o';
            while ischar(tline)
                counter = counter + 1;
                tline = fgets(fid);
                if tline ~= -1
                    B{counter} = strsplit(tline, ' ', 'CollapseDelimiters',true);
                else
                end
            end
            m=0;
            for i=1:length(B)
                if strcmp(Element.Name, B{i}{1})
                    m = i;
                else
                end
            end   
            sellmeier = B{m}{2};
            A(l) = evalin(symengine, sellmeier);
            fclose(fid);
            c=2.9979e8;
            w = sym('w');
            S(w) = A(((c*2*pi)/w));
            S1 = diff(S);
            GD_s=Element.Length*(S+w*(S1))/c;
            GD=matlabFunction(GD_s);
            S2 = diff(S,2);
            GDD_s = (Element.Length/c)*(2*S1+w*S2);
            GDD = matlabFunction(GDD_s);
            coeff1 = GD((c*2*pi)./Element.CentralWavelength);
            coeff2 = GDD((c*2*pi)./Element.CentralWavelength);
            dispersionCoeff=[coeff1,coeff2];
        end
     end
end    