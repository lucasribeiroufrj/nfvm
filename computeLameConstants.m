function [ lam1, lam2 ] = computeLameConstants( params )
%computeLameConstants 
%   Given a struct like this struct('E', 1e+7, 'Nu', 0.3), it returns
% the Lame's constants Lambda and Mu.
%   Given a struct like this struct('Lambda', 1e+11, Mu', 1e+7), 
% it returns the Lame's constants E and Nu.

if isfield(params, 'Lambda') && isfield(params, 'Mu')
    
    Lambda = params.Lambda;
    Mu = params.Mu;
    
    E = Mu * (3*Lambda + 2*Mu) / (Lambda + Mu);
    Nu = Lambda / (2*(Lambda + Mu));
    
    lam1 = E;
    lam2 = Nu;
    
elseif isfield(params, 'E') && isfield(params, 'Nu')
    
    E = params.E;
    Nu = params.Nu;
    isPlaneStress = params.planeStress;

    if isPlaneStress
        Lambda = E*Nu / (1-Nu*Nu);
    else
        Lambda = E*Nu / ( (1+Nu)*(1-2*Nu) );
    end

    Mu = E / (2*(1+Nu));
    
    lam1 = Lambda;
    lam2 = Mu;
    
else
    
    ss = [ 'struct(''E'', 1e+7, ''Nu'', 0.3)', ' or ', ...
        'struct(''Lambda'', 1e+11, Mu'', 1e+7)'];
    
    ME = MException('MyComponent:inputError',...
                'Invalid input. type something like: \n%s', ss);
    throw(ME)
    
end

end

