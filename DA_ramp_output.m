function output = DA_ramp_output(a,e,mI,da,type,gain,varargin)

%DA_RAMP_OUTPUT Dopamine affected unit output 
%   O = DA_RAMP_OUTPUT(A,E,M,DA,T,G,P) computes the output O of a LI unit layer,
%   given activation A, threshold E, and initial slope M arrays, dopamine level DA,
%   dopamine receptor type T (1 = D1, 2 = D2), and gain G. The pivot point
%   P is a optional parameter in that it only needs to be specified for D1
%   receptors.
%
%   Reference: Humphries, M.D. (2003). High level modeling of dopamine mechanisms 
%   in striatal neurons. Technical Report ABRG 3. Dept. Psychology
%   University of Sheffield, UK.
% 
%   Mark Humphries 21/1/2005

% check for pivot
if type==1 & nargin < 7
    error('Must specify pivot parameter for D1 receptors')
end

%%%%%%%%%%%%%
%%% below is an optimised verison of this for arrays of activations...
% if a < e
%     output = 0;
% elseif a <= 1/m + e
%     output = m * (a - e);
% else 
%     output = 1;
% end

if type == 1    % D1 model
    p = varargin{1};
    m = mI + gain .* da; 
    
    % classify outputs
    limit = (1 - (1 - m) .* p) ./ m + e;
	case_zero = find(a < e);
	case_a = find(a >= e & a <= limit);
	case_one = find(a > limit);

	% compute outputs
	output(case_zero) = 0;
	output(case_a) = m .* (a(case_a) - e) + (1 - m) .* p;
	output(case_one) = 1;
    
elseif type == 2 % D2 model
    m = mI - gain .* da; 
    
    % case statement same as original ramp function
	case_zero = find(a < e);
	case_a = find(a >= e & a <= 1/m +e);
	case_one = find(a > 1/m + e);
	
	output(case_zero) = 0;
	output(case_a) = m .* (a(case_a) - e);
	output(case_one) = 1;
else
    error('Unknown dopamine receptor type specified')
end

    

