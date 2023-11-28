MATLAB version of the extended basal ganglia model

Files:
extended_model_ZOH.m
HG_engine.m
ramp_output.m
DA_ramp_output.m


Notes:
(1) extended_model_ZOH.m is a MATLAB script that replicates the Simulink version

(2) HG_engine.m	is a MATLAB function that encapsulates multiple versions of the model: read its Help comments for full details.
 
(3) ramp_output.m and DA_ramp_output.m are two different forms of the piece-wise linear output function. The latter is a modified form given in (Humphries, 2003) that captures the effects of dopamine on striatal neuron output. 


References:
Humphries, M. D. & Gurney, K. N. (2002). The role of intra-thalamic and thalamocortical circuits in action selection. Network: Computation in Neural Systems, 13, 131-156.

Humphries, M.D. (2003). High level modeling of dopamine mechanisms in striatal neurons. Technical Report ABRG 3. Dept. Psychology University of Sheffield, UK.
 