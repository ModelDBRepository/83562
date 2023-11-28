function [winner,A,O,step_counter] = HG_engine(saliences,DA_sel,DA_cont,dt,tolerance,max_steps,theta,varargin)

% HG_ENGINE solution engine for the Humphries & Gurney (2002) extended BG model 
%
%	HG_ENGINE(S,DA1,DA2,DT,TOLERANCE,MAX_STEPS,THETA) where S is an array of salience values for the actions
%       represented on the BG channels (the number of channels is implicitly defined by the length of the 
%       salience array). Dopamine levels in the selection and control pathway are set by the values of DA1 and DA2.
%       The time-step is DT, and the model is run until either the change in activation of all units is less than TOLERANCE
%       or MAX_STEPS has been reached. If the output of any GPi channel is below THETA then that channel is considered
%       selected. Returns: an array of the winning action(s) (channel(s)), or the empty matrix [] if no winner
%
%	HG_ENGINE(...,SWITCH) where SWITCH = 'hard' enforces hard switching so that a maximum of one selection is 
%       made. When more than one GPi unit's output is below THETA then the channel of the lowest is returned, else [] is returned.
%       Where SWITCH = 'gate', returns a vector containing the proportion
%       of output below THETA for each channel, where 0 indicates that the
%       output is above THETA, and 1 indicates no output.
%
%   [W,nA,nO,STEPS] = HG_ENGINE(...,A,O) are the matrices of activations A and outputs O of
%   all the units from the previous competition. By column: [SD1 SD2 STN
%   GPe GPi]. W is the array of winner(s) if any. Specifying nA and nO returns the corresponding matrices from
%   the current simulation to re-use as arguments for the next one. Put
%   SWITCH = [] if no need to specify this parameter. Will also return
%   STEPS, the number of steps to convergence.
%
%   HG_ENGINE(...,FLAG) any combination of the following options creates (set A=[], O=[] if not required):
%       'g': includes the connections from the Gurney et al. (2004) Network
%       paper too
%
%       'd': includes the new DA model from my technical report: Humphries, M.D. (2003). High level %	     modeling of dopamine mechanisms in striatal neurons. ABRG 3. Dept. Psychology University of Sheffield, UK.
%
%  See notes in GPR_ENGINE. In addition: The model implemented here contains two modifications from the published model
%       (1) There is no GPi->TRN connection (due to lack of evidence)
%       (2) The STN connection weights are taken from the GPR model
%
%   Author: Mark Humphries 21/1/05

%%% MODEL PARAMETERS
NUM_CHANNELS = length(saliences);
NUM_NUCLEI = 8;

%% weight values as defined by HG
W_Sal_SEL = 0.5;
W_Sal_CONT = 0.5;
W_Sal_STN = 0.5;
W_Sal_MCtx = 1;

W_MCtx_SEL = 0.5;
W_MCtx_CONT = 0.5;
W_MCtx_STN  = 0.5;
W_MCtx_VL = 1;
W_MCtx_TRN = 1;

W_VL_MCtx = 1;
W_VL_TRN = 1;

W_TRN_between = -0.7;
W_TRN_within = -0.1;

W_SEL_GPi = -1;
W_CONT_GPe = -1;
W_STN_GPi = 0.9;
W_STN_GPe = 0.9;
W_GPe_STN = -1;
W_GPe_GPi = -0.3;
W_GPi_VL = -1;

%% additional weights are zero by default
W_GPi_GPi = 0;
W_GPe_GPe = 0;
W_SEL_GPe = 0;

%% thresholds as defined by GPR / HG
e_SEL = 0.2;
e_CONT = 0.2;
e_STN = -0.25;
e_GPe = -0.2;
e_GPi = -0.2;
e_MCtx = 0;
e_VL = 0;
e_TRN = 0;

%%% INITIALISE ARRAYS

% activity arrays
A = zeros(NUM_CHANNELS,NUM_NUCLEI); % 1 = MCtx, 2 = VL, 3 = TRN, 4 = SD1, 5 = SD2, 6 = STN, 7 = GPe, 8 = GPi
old_A = zeros(NUM_CHANNELS,NUM_NUCLEI);
delta_a = ones(NUM_CHANNELS,NUM_NUCLEI);

% output arrays
O = zeros(NUM_CHANNELS,NUM_NUCLEI);

%%% Optional parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type = 'soft';
blnNewDA = 0;
gain_DA_sel = DA_sel;      
gain_DA_cont = DA_cont;

if nargin >= 8 & ~isempty(varargin{1}) type = varargin{1}; end
if nargin >= 9 & ~isempty(varargin{2})
    [rA cA] = size(varargin{2});    
    if cA ~= NUM_NUCLEI error(['Activation matrix must have' num2str(NUM_NUCLEI) ' columns for HG model']); end
    if rA ~= NUM_CHANNELS error('Activation matrix must have the same number of rows as specified saliences'); end
    A = varargin{2};
end
if nargin >= 10 & ~isempty(varargin{3})
    [rO cO] = size(varargin{3});    
    if cO ~= NUM_NUCLEI error(['Output matrix must have' num2str(NUM_NUCLEI) 'columns for HG model']); end
    if rO ~= NUM_CHANNELS error('Output matrix must have the same number of rows as specified saliences'); end
    O = varargin{3};
end
if nargin >= 11 ~isempty(varargin{4})
    if findstr(varargin{4},'g') % include weights from Gurney et al (2004) Network paper
        W_GPi_GPi = -0.2;
		W_GPe_GPe = -0.2;
		W_SEL_GPe = -0.25;   
    end
    if findstr(varargin{4},'d') % include new DA model 
        gain_DA_sel = 0;        % set gain DA to zero
        gain_DA_cont = 0;
        blnNewDA = 1;       % set flag to use new output functions
        % parameter values from tech. report
        e_SEL = 0.1;
        e_CONT = 0.1;       
        gain_SEL = 0.8;
        gain_CONT = 0.8;
        pivot = 0.1;
    end
end

%%% ARTIFICAL UNIT PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 25;                     % gain
m = 1;                      % slope

decay_constant = exp(-k*dt);    
int_vec = ones(NUM_CHANNELS,1);

%%% SIMULATE MODEL
step_counter = 0;

[row col] = size(saliences);
if row < col
    c = saliences';
else
    c = saliences;
end


while step_counter < max_steps & sum(sum(delta_a > tolerance)) > 0 
    step_counter = step_counter + 1;
    old_A = A;
    
    %% MOTOR CORTEX
    u_MCtx = c .* W_Sal_MCtx + O(:,2) .* W_VL_MCtx;
    A(:,1) = (A(:,1) - u_MCtx) * decay_constant + u_MCtx;
    
    %% VL thalamus
    temp = (sum(O(:,3)) .* int_vec) - O(:,3);     %% removes own input from each summed value
    u_VL = O(:,1) .* W_MCtx_VL + O(:,2) .* W_GPi_VL + O(:,3) .* W_TRN_within + temp .* W_TRN_between;
    A(:,2) = (A(:,2) - u_VL) * decay_constant + u_VL;

    %% TRN
    u_TRN = O(:,1) .* W_MCtx_TRN + O(:,2) .* W_VL_TRN;
    A(:,3) = (A(:,3) - u_TRN) * decay_constant + u_TRN;
    
    %% STRIATUM D1
    u_SEL = (c .* W_Sal_SEL + O(:,1) .* W_MCtx_SEL) .* (1 + gain_DA_sel);
    A(:,4) = (A(:,4) - u_SEL) * decay_constant + u_SEL;
    
    %% STRIATUM D2
    u_CONT = (c .* W_Sal_CONT + O(:,1) .* W_MCtx_CONT) .* (1 - gain_DA_cont);
    A(:,5) = (A(:,5) - u_CONT) * decay_constant + u_CONT;

    %% STN
    u_STN = c .* W_Sal_STN + O(:,1) .* W_MCtx_STN + O(:,7) .* W_GPe_STN;
    A(:,6) = (A(:,6) - u_STN) * decay_constant + u_STN;
    
    %% GPe
    temp = (sum(O(:,7)) .* int_vec) - O(:,7);     %% removes own input from each summed value
    u_GPe = sum(O(:,6)) .* W_STN_GPe + O(:,5) .* W_CONT_GPe + O(:,4) .* W_SEL_GPe + temp .* W_GPe_GPe;;
    A(:,7) = (A(:,7) - u_GPe) * decay_constant + u_GPe;

    %% GPi
    temp = (sum(O(:,8)) .* int_vec) - O(:,8);     %% removes own input from each summed value    
    u_GPi = sum(O(:,6)) .* W_STN_GPi + O(:,7) .* W_GPe_GPi + O(:,4) .* W_SEL_GPi + temp .* W_GPi_GPi;
    A(:,8) = (A(:,8) - u_GPi) * decay_constant + u_GPi;
    
    % calculate outputs
    O(:,1) = ramp_output(A(:,1),e_MCtx,m)';
    O(:,2) = ramp_output(A(:,2),e_VL,m)'; 
    O(:,3) = ramp_output(A(:,3),e_TRN,m)'; 
    if blnNewDA
        O(:,4) = DA_ramp_output(A(:,4),e_SEL,m,DA_sel,1,gain_SEL,pivot)';    
        O(:,5) = DA_ramp_output(A(:,5),e_CONT,m,DA_cont,2,gain_CONT)';    
    else
        O(:,4) = ramp_output(A(:,4),e_SEL,m)';    
        O(:,5) = ramp_output(A(:,5),e_CONT,m)';    
    end
    O(:,6) = ramp_output(A(:,6),e_STN,m)'; 
    O(:,7) = ramp_output(A(:,7),e_GPe,m)';     
    O(:,8) = ramp_output(A(:,8),e_GPi,m)';
    
    % change in activation for each unit
    delta_a = abs(A - old_A);
end

% determine winner(s)
winner = [];
switch type
case 'soft'
    winner = find(O(:,8) < theta);
case 'hard'
    temp = find(O(:,8) < theta);
    
    if ~isempty(temp) winner = find(O(:,8) == min(O(temp,8))); end
    
    if length(winner) > 1
        winner = [];    % can only return one winner
    end
case 'gate'
    winner = zeros(NUM_CHANNELS,1);
    output = O(:,8); 
    idxs = find(output < theta);
    winner(idxs) = 1 - (theta - output(idxs));    
end
