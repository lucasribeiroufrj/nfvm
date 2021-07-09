classdef SteadyTime < Time & handle
    %SteadyTime Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        hasAlreadyRun_ = false;
    end
    
    methods
        
        function obj = SteadyTime()
            
            obj@Time(0, 0, 0, 0);
        end
        
        function obj = increment(obj)
            %increment Do nothing.
            
        end
        
        function obj = pushDeltaT(obj, deltaT)
            
            error('Invalid operation');
        end
        
        function obj = changeDelta_t(obj, deltaT)
            
            error('Invalid operation');
        end
        
        function shouldRun = run(obj)
            %run Return true if run should continue.
            
            if ~obj.hasAlreadyRun_;
                
                shouldRun = true;
                obj.hasAlreadyRun_ = true;
            else
                shouldRun = false;
            end
        end
    end
    
end

