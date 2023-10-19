function Dnext = forward_policy(D, a_i, a_pi,Pi)
%--------------------------------------------
% Description: function to do forward iteration of stationary
% distribution D(e,a)
%--------------------------------------------
% Input
% Dt:   Current stationary distribution D(e,a)
% a_i:  Lottery values of policy function a'(e,a)
% a_pi: Lottery probability of asset
% Pi:   Transition matrix of income
%----------------------------------------------
% Output
% Dt+1: Next period stationary distribution Dt+1(e',a')
%----------------------------------------------

% Preallocate the resulting matrix
    [m, n] = size(D);
    Dend = zeros(m, n);
    for e = 1:m
        for a = 1:n
            % Send pi(e,a) of the mass to gridpoint i(e,a)
            Dend(e, a_i(e,a)+1) = Dend(e, a_i(e,a)+1) + a_pi(e,a)*D(e,a);
            % Check bounds to ensure we don't exceed matrix dimensions
            if a_i(e, a) + 2 <= n
                % Send 1-pi(e,a) of the mass to gridpoint i(e,a)+1
                Dend(e, a_i(e,a)+2) = Dend(e, a_i(e,a)+2) + (1-a_pi(e,a))*D(e,a);
            end
        end
    end
    % Dend is the "end-of-period distribution" D_end(e,a');
    % Dt+1(e',a')=Pi'*D_end(e,a');
    Dnext=Pi' * Dend;
end