%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  DATE: 20/10/2003
%%%%  WHAT: M-code (script) version of the extended model
%%%%		Parameters are from Humphries & Gurney (2002) Network, 13, 131-156.
%%%%  AUTHOR: Mark Humphries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%%% MODEL PARAMETERS
NUM_CHANNELS = 6;
da_sel = 0.2;     % dopamine level
da_cont = 0.2;

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
W_STN_GPi = 0.8;
W_STN_GPe = 0.8;
W_GPe_STN = -1;
W_GPe_GPi = -0.4;
W_GPi_VL = -1;

e_SEL = 0.2;
e_CONT = 0.2;
e_STN = -0.25;
e_GPe = -0.2;
e_GPi = -0.2;
e_MCtx = 0;
e_VL = 0;
e_TRN = 0;

%%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 6;                      % length of simulation
dt = 0.001;                  % time-step
sim_steps = t / dt;

% activity arrays
a_MCtx = zeros(NUM_CHANNELS,1);
a_VL = zeros(NUM_CHANNELS,1);
a_TRN = zeros(NUM_CHANNELS,1);
a_SEL = zeros(NUM_CHANNELS,1);
a_CONT = zeros(NUM_CHANNELS,1);
a_STN = zeros(NUM_CHANNELS,1);
a_GPe = zeros(NUM_CHANNELS,1);
a_GPi = zeros(NUM_CHANNELS,1);

% output arrays
o_MCtx = zeros(NUM_CHANNELS,sim_steps);
o_VL = zeros(NUM_CHANNELS,sim_steps);
o_TRN = zeros(NUM_CHANNELS,sim_steps);
o_SEL = zeros(NUM_CHANNELS,sim_steps);
o_CONT = zeros(NUM_CHANNELS,sim_steps);
o_STN = zeros(NUM_CHANNELS,sim_steps);
o_GPe = zeros(NUM_CHANNELS,sim_steps);
o_GPi = zeros(NUM_CHANNELS,sim_steps);

%%% SALIENCE INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch1_onset = 1;
ch1_size = 0.4;

ch2_onset = 3;
ch2_size = 0.6;

transient_onset = 4;
transient_off = 6;
transient_size = 0.2;

c = zeros(NUM_CHANNELS,1);

%%% ARTIFICAL UNIT PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 25;                     % gain
m = 1;                      % slope

decay_constant = exp(-k*dt);    
trn_vec = ones(NUM_CHANNELS,1);

tic
%%% SIMULATE MODEL
for steps = 2:sim_steps
    
    %% calculate salience changes
    if steps * dt == ch1_onset
        c(1) = ch1_size;
    end
    if steps * dt == ch2_onset
        c(2) = ch2_size;
    end
    if steps * dt == transient_onset
        c(1) = c(1) + transient_size;
    end
    if steps * dt == transient_off
        c(1) = ch1_size;
    end
    
    %%% OUTPUTS - calculated from previous time-step
    o_MCtx(:,steps) = ramp_output(a_MCtx,e_MCtx,m)';    
    o_VL(:,steps) = ramp_output(a_VL,e_VL,m)'; 
    o_TRN(:,steps) = ramp_output(a_TRN,e_TRN,m)'; 
    o_SEL(:,steps) = ramp_output(a_SEL,e_SEL,m)';  
    o_CONT(:,steps) = ramp_output(a_CONT,e_CONT,m)';    
    o_STN(:,steps) = ramp_output(a_STN,e_STN,m)'; 
    o_GPe(:,steps) = ramp_output(a_GPe,e_GPe,m)';     
    o_GPi(:,steps) = ramp_output(a_GPi,e_GPi,m)';

    %% Motor Cortex
    u_MCtx = c .* W_Sal_MCtx + o_VL(:,steps) .* W_VL_MCtx;
    a_MCtx = (a_MCtx - u_MCtx) * decay_constant + u_MCtx;
    
    %% VL thalamus
    temp = (sum(o_TRN(:,steps)) .* trn_vec) - o_TRN(:,steps);     %% removes own input from each summed value
    u_VL = o_MCtx(:,steps) .* W_MCtx_VL + o_GPi(:,steps) .* W_GPi_VL + o_TRN(:,steps) .* W_TRN_within + temp .* W_TRN_between;
    a_VL = (a_VL - u_VL) * decay_constant + u_VL;
    
    %% TRN
    u_TRN = o_MCtx(:,steps) .* W_MCtx_TRN + o_VL(:,steps) .* W_VL_TRN;
    a_TRN = (a_TRN - u_TRN) * decay_constant + u_TRN;

    %% STRIATUM D1
    u_SEL = (c .* W_Sal_SEL + o_MCtx(:,steps) .* W_MCtx_SEL) .* (1 + da_sel);
    a_SEL = (a_SEL - u_SEL) * decay_constant + u_SEL;

    %% STRIATUM D2
    u_CONT = (c .* W_Sal_CONT + o_MCtx(:,steps) .* W_MCtx_CONT) .* (1 - da_cont);
    a_CONT = (a_CONT - u_CONT) * decay_constant + u_CONT;

    %% STN
    u_STN = c .* W_Sal_STN + o_MCtx(:,steps) .* W_MCtx_STN + o_GPe(:,steps) .* W_GPe_STN;
    a_STN = (a_STN - u_STN) * decay_constant + u_STN;
    
    %% GPe
    u_GPe = sum(o_STN(:,steps)) .* W_STN_GPe + o_CONT(:,steps) .* W_CONT_GPe;
    a_GPe = (a_GPe - u_GPe) * decay_constant + u_GPe;

    %% GPi
    u_GPi = sum(o_STN(:,steps)) .* W_STN_GPi + o_GPe(:,steps) .* W_GPe_GPi + o_SEL(:,steps) .* W_SEL_GPi;
    a_GPi = (a_GPi - u_GPi) * decay_constant + u_GPi;
    
end

toc

figure(1)
clf
plot(o_GPi(1,:),'r')
hold on
plot(o_GPi(2,:),'b')
