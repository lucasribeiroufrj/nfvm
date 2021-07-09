classdef OscillationMonitor < handle
    %OscillationMonitor Detects when solution stops changing
    %   Detailed explanation goes here
    
    properties
        
        repetition_;
        maxNumberOfRepetitions_ = 15;
        repArray_;
        threshold_ = 1e-2;
    end
    
    methods
        function obj = OscillationMonitor(maxNumberOfRepetitions, threshold)
            
            if nargin == 2
                obj.maxNumberOfRepetitions_ = maxNumberOfRepetitions;
                obj.threshold_ = threshold;
            end
            
            obj.restart();
        end
        
        function obj = restart(obj)
            
            obj.repetition_ = 0;
            obj.repArray_ = zeros(1, obj.maxNumberOfRepetitions_);
        end
        
        function shouldStop = stop(obj, residual)
            
            shouldStop = false;
            
            rep = mod(obj.repetition_, obj.maxNumberOfRepetitions_);

            if obj.repetition_ + 1 >= obj.maxNumberOfRepetitions_

                obj.repArray_(rep+1) = residual;

                res = sum(obj.repArray_) / obj.maxNumberOfRepetitions_;

                if residual ~= 0 ...
                    && (abs(res-residual) / abs(residual)) < obj.threshold_

                    obj.printAlert();
                    shouldStop = true;
                    
                    return;
                end
            else
                obj.repArray_(rep+1) = residual;
            end
            
            obj.repetition_ = obj.repetition_ + 1;
        end
        
        function obj = printAlert(obj)
            n = obj.maxNumberOfRepetitions_;
    disp('┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓')
    disp('┃    ATTENTION: Solution has not been changing for    ┃');
 fprintf('┃               %2d step(s).                           ┃\n', n);
    disp('┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛')
        end
    end
    
end

