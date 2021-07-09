classdef Time < handle
    %Time Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        timeIndex_
        startTime_
        endTime_
        deltaT_
        t_
        t_o
        deltaT_o_
        queueDeltaT_ = []
    end
    
    methods
        
        function obj = ...
            Time...
            (...
                startTime,...
                endTime,...
                deltaT,...
                deltaT_o...
            )
            
            obj.timeIndex_ = 0;
            obj.startTime_ = startTime;
            obj.endTime_ = endTime;
            obj.deltaT_ = deltaT;
            obj.t_ = obj.startTime_;
            obj.deltaT_o_ = deltaT_o;
            obj.t_o = obj.t_ - obj.deltaT_o_;
        end
        
        %% BEG - Data member access
        function timeIndex = timeIndex(obj)
            
            timeIndex = obj.timeIndex_;
        end
        
        function startTime = startTime(obj)
           
            startTime = obj.startTime_;
        end
        
        function endTime = endTime(obj)
           
            endTime = obj.endTime_;
        end
        
        function obj = setEndTime(obj, newEndTime)
            
            obj.endTime_ = newEndTime;
        end
        
        function deltaT = deltaT(obj)
            
            deltaT = obj.deltaT_;
        end
        
        function time = value(obj)

            time = obj.t_;
        end
        
        function time = valueOld(obj)

            time = obj.t_o;
        end
        
        function time = valueOldOld(obj)

            time = obj.t_o - obj.deltaT_o_;
        end
        
        function bool = isFirstTimeStep(obj)
           bool = obj.timeIndex_ == 1; 
        end
        %% END
        
        function obj = increment(obj)
            %increment Go to the next time step.
            
            obj.deltaT_o_ = obj.deltaT_;
            
            if ~isempty(obj.queueDeltaT_)
                
                obj.deltaT_ = obj.queueDeltaT_(end);
                obj.queueDeltaT_ = obj.queueDeltaT_(1:end-1);
            end
            
            obj.t_o = obj.t_;
            obj.t_ = obj.t_ + obj.deltaT_;
            obj.timeIndex_ = obj.timeIndex_ + 1;
        end
        
        function obj = pushDeltaT(obj, deltaT)
            
            obj.queueDeltaT_ = [ deltaT obj.queueDeltaT_ ];
        end
        
        function obj = changeDelta_t(obj, deltaT)
            
            obj.deltaT_ = deltaT;
        end
        
        function obj = resetTimeIndex(obj)
            
            obj.timeIndex_ = 0;
        end
        
        function value = run(obj)
            %run Return true if run should continue.
            
            value =  obj.t_ < (obj.endTime_ - 0.5*obj.deltaT_);
        end
        
        function value = deltaT_o(obj)
           
            value = obj.deltaT_o_;
        end
    end
    
end

