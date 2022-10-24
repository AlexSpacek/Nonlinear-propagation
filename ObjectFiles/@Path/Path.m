classdef Path < handle  
    properties
        Name
    end   
    methods (Abstract) %all of the subclasses must have the method defined (can be defined differently)
        dispersionCoeff = getDispersionCoefficients(thisField)
    end    
    methods
        function thisPath=Path(name)  %all the shapes share this exact function
            if nargin==1
                if ischar(name)
                    thisPath.Name=name;
                else
                    error('Name of the element is not a string')
                end
%                 if isnumeric(order)
%                     thisPath.Order=order;
%                 else
%                     error('Order of the element must be an integer value')
%                 end   
            else 
                error('Must provide name of the material used');
            end
        end
    end
end