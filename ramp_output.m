%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DATE: 7/7/2003
%%%% WHAT: piece-wise linear output function (i.e. ramp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = ramp(a,e,m)

%%%%%%%%%%%%%
%%% below is an optimised verison of this for arrays of activations...
% if a < e
%     output = 0;
% elseif a <= 1/m + e
%     output = m * (a - e);
% else 
%     output = 1;
% end

case_zero = find(a < e);
case_a = find(a >= e & a <= 1/m +e);
case_one = find(a > 1/m + e);

output(case_zero) = 0;
output(case_a) = m .* (a(case_a) - e);
output(case_one) = 1;

