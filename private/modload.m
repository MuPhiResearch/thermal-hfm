function modload
%Initialize Thermal EDFM Module
% Couples the hfm (EDFM) module with the geothermal module for thermal
% simulation through fractured reservoirs.
%
% Dependencies:
%   hfm           - EDFM grid/NNC preprocessing (2D and 3D pipelines)
%   geothermal    - thermal transport models and state functions
%   compositional - provides CompositionalMixture, the superclass of
%                   CompositionalBrineFluid. GeothermalModel (parent of
%                   GeothermalHFMModel) unconditionally constructs a
%                   CompositionalBrineFluid in its constructor, so
%                   compositional must be on the path for the model to
%                   instantiate.
%
% Optional (loaded by the relevant examples, not here, to keep EDFM-only
% usage lightweight):
%   shale - required for pEDFM (fractureMethod='pedfm') and the pEDFM preprocessing
%           (setupPEDFMOpsTPFA, transmultpEDFM, pMatFracNNCs3D, EDFMshalegrid).
%   dfm   - required by the DFM-vs-EDFM comparison example.

mrstModule add hfm geothermal compositional;

end
