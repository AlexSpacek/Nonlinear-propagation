classdef Lens < Path
    properties
        Length
        NonlinearIndex
        RefractiveIndex
        RefractiveIndexVector
        WavelengthVector
        Dz
        CentralWavelength
        RungeKuttaSteps
        Focus
    end


    methods
        function thisLens=Lens(name,length2,numberOfPoints,rungeKuttaSteps,nonlinearIndex,wavelengthVector,centralWavelength,focus) 
            thisLens@Path(name)
            if nargin==8 
                if isnumeric(length2)
                    thisLens.Length=length2;
                else
                    error('Length must be numeric')
                end
                if isnumeric(rungeKuttaSteps)
                    thisLens.RungeKuttaSteps=rungeKuttaSteps;
                else
                    error('Number of steps for Runge Kutta must be numeric')
                end   
                if isnumeric(nonlinearIndex)
                    thisLens.NonlinearIndex=nonlinearIndex;
                else
                    error('Nonlinear index must be numeric')
                end  
                if isnumeric(wavelengthVector)
                    thisLens.WavelengthVector=wavelengthVector;
                else
                    error('Wavelength vector must be numeric')
                end
                if isnumeric(centralWavelength)
                    thisLens.CentralWavelength=centralWavelength;
                else
                    error('Central wavelength must be numeric')
                end 
                if isnumeric(numberOfPoints)
                    thisLens.Dz=thisLens.Length/numberOfPoints;
                else
                    error('The number of sampling points must be numeric')
                end
                if isnumeric(focus)
                    thisLens.Focus=focus;
                else
                    error('Focal distance must be numeric')
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
                thisLens.RefractiveIndexVector = double(A(wavelengthVector));
                thisLens.RefractiveIndex = double(A(centralWavelength));
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
            c = 3e8;
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