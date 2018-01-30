classdef myPoint
    %MYPOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess='public', SetAccess='public')
    x;
    y;
    end
    
    methods 
        function obj = myPoint(x,y)
           obj.x = x;
           obj.y = y;
        end
    end
    
end

