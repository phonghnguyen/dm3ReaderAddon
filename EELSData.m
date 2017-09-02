classdef EELSData
    properties
        %writen on import
        Data
            %properties of the dimentions of 'Data'
            dim 
            origin
            scale
            unit
        expo
        apat
        disp
        mag
        %dependent proporties
        %writen by later functions
        shift
    end
    methods
        function dim = get.dim(obj)
            %generates a vector contaning the dimentions of the data
            dim = size(obj.Data);
        end
        function out = axis(obj,input)
            %Given a dimention in 'input' generates a row vector containing
            %the indepentent vatiables of 'Data' along that axis. i.e.
            %obj.axis(3) is a row vector containing the energy values of
            %the EELS data in 'Data'
            out = ((0:obj.dim(input)-1)+obj.origin(input)).*obj.scale(input);
        end
        function out = saxis(obj,xy)
            %Given an input 'xy' of the form "xy = [x,y]" where x and y are
            %indicies of elements in the first 2 dimentions of 'Data'
            %generates a vector containing the 
            if isempty(obj.shift)
                error('Must first find zero loss shift data')
            else
            out = obj.axis(3) + obj.shift(xy(1),xy(2));
            end
        end
    end
end
