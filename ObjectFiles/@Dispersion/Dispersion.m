classdef Dispersion < Path
    properties
        GDD
        TOD
        FOD
    end


    methods
        function thisDispersion=Dispersion(name,GDD,TOD,FOD)
            thisDispersion@Path(name)
            if nargin==4 
                if isnumeric(GDD)
                    thisDispersion.GDD=GDD;
                else
                    error('GDD must be numeric')
                end
                if isnumeric(TOD)
                    thisDispersion.TOD=TOD;
                else
                    error('TOD must be numeric')
                end
                if isnumeric(FOD)
                    thisDispersion.FOD=FOD;
                else
                    error('FOD must be numeric')
                end
            end
        end
       
     function  dispersionCoeff = getDispersionCoefficients(Element) 
         dispersionCoeff=[Element.GDD,Element.TOD,Element.FOD];
     end
    end  
end