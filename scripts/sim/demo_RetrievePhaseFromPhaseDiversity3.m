%%
% File: demo_RetrievePhaseFromPhaseDiversity3.m
% Author: Santiago Echeverri, Ren� Restrepo, Carlos Cuartas and N�stor
% Uribe
% Date:	  19/07/2016
% Modification:	  10/07/2020
% Vortex-enhanced coherent-illumination phase diversity for phase retrieval in coherent imaging systems 
% Optics Letters Vol. 41, Issue 8, pp. 1817-1820 (2016). https://doi.org/10.1364/OL.41.001817
%%
%  This program retrieves phase by using the proposed coherent phase diversity.
%%
%
%		This file contains one example of the PD implementation in this case,
%		the procedure is employed using four diversities, two with topological
%		charge [-2,2] and no Zernikes and two with topological charge [0, 0] and
%		[-0.5, 0.5] lambda astigmatism
%
%		NOTE: This program does ignore piston and tilts.
%
%%

%% Packages
% Used packages (for octave purposes)
% pkg load image;
% pkg load parallel;
% pkg load optim;
% warning("off");

%% Constants
FTSize = 1024;% Images size
gridSize = 120;% Padarray size
gaussianC = 0.3;% Gaussian beam diameter, this is with respect to gridsize, in this case gaussianC = 0.3 * gridSize
rLimit = 74/ gridSize;% FT gaussian effective size
verbosity = 1;
dx = 0;
dy = 0;
%mType = 'nelder';
mType = 'nonlin';
mPars = {1000 1e-12 4000 100};% Minimization parameters, for this procedure [isz, tolFun]
% The isz is the initial step size which determines the magnitude of the steps in first iteration - recommended 1-
% tolFun is the tolerance of the function determined by the change between two successive iteration, if the error change is
% smaller than the given tolerance, the algorithm stops. -recommended range [1e-4 - 1e-7]-.
addpath(genpath('../../matlab'));


%% Diversities
% We want four diversities, so four lDiv are declarated, as mentioned initially, those are [0, 0, -2, 2]. The Zernike $kDiv$ are
% [0 0.5], [0 -0.5], [0 0], [0 0], as implementated
kDiv = {[0 0.5 0 0 0 0 0 0 0 0 0 0 0]', [0 -0.5 0 0 0 0 0 0 0 0 0 0 0]', [0 0 0 0 0 0 0 0 0 0 0 0 0]', [0 0 0 0 0 0 0 0 0 0 0 0 0]'};
lDiv = [0, 0, -2, 2];

%% Desired aberration
% In our notation, aberrations follow
% [defocus, astigmatism, astigmatism45°, comaX, comaY, spherical, trefoil, trefoil45°, ...
% secondAstigmatism, secondAstigmatism45°, secondComaX, secondComaY, ...],
ZAberration = [0 0.3 0 -0.4 0 0 0.5 0 0 0 0 0 0]'; % This is astigmatism

%% Image range
colRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;
rowRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;

%% Diversity images creation
% We use CalcOAMBeamFTFromAberrations2 which is the image simulator in order to produce the desire diversities.
% This one has OAM 0, 0.5ast
coeff = 3:size(ZAberration,1)+2;
order = lDiv(1);
ZBaseAberration = kDiv{1} + ZAberration;
[~, IntensityCenterZ3( :, :, 1), ~, ~, ~, ~, ~] = CalcOAMBeamFTFromAberrations2(FTSize, gridSize, gaussianC, rLimit, ZBaseAberration, dx, dy, order, coeff);

% This one has OAM 0, -0.5ast
order = lDiv(2);
ZBaseAberration = kDiv{2} + ZAberration;
[~, IntensityCenterZ3( :, :, 2), ~, ~, ~, ~, ~] = CalcOAMBeamFTFromAberrations2(FTSize, gridSize, gaussianC, rLimit, ZBaseAberration, dx, dy, order, coeff);

% This one has OAM -2, 0ast
order = lDiv(3);
ZBaseAberration = kDiv{3} + ZAberration;
[~, IntensityCenterZ3( :, :, 3), ~, ~, ~, ~, ~] = CalcOAMBeamFTFromAberrations2(FTSize, gridSize, gaussianC, rLimit, ZBaseAberration, dx, dy, order, coeff);

% This one has OAM 2, -0.5ast
order = lDiv(4);
ZBaseAberration = kDiv{4} + ZAberration;
[~, IntensityCenterZ3( :, :, 4), ~, ~, ~, ~, ~] = CalcOAMBeamFTFromAberrations2(FTSize, gridSize, gaussianC, rLimit, ZBaseAberration, dx, dy, order, coeff);

% Only keep the area with information
ROIImageZ3 = IntensityCenterZ3( colRangeSLM, rowRangeSLM, :);
clear IntensityCenterZ3

%% PD run
% For experimental images replace ROIImageZ3 with your images organized as the diversities

% Initial seed for the algorithm
ZBaseAberration = zeros(15,1);
coeff = 3:size(ZBaseAberration,1)+2;

tic; % time
[ZAberrationsZ3, ~, ~, ~, ~, ~, ~, ~, waveAberationPD, PDIntensityZ3] = RetrievePhaseFromPhaseDiversity1(ROIImageZ3, FTSize, gridSize, gaussianC, rLimit, ZBaseAberration, 0, dx, dy, kDiv, lDiv, mType, mPars, verbosity, coeff); % Run phase diversities

% note that printed in terminal is the final error and the aberrations retrieved, which clearly is 1lambda astigmatism
tiempoPDZ3 = toc;

% Final result info
errPDZ3 = sum(sum(abs(ROIImageZ3(:,:,1) - PDIntensityZ3(colRangeSLM,rowRangeSLM,2)) / (size(PDIntensityZ3,1) * size(PDIntensityZ3,2)), 1), 2);
disp(['Error final por PDZ3 ',num2str(errPDZ3),' tiempo ', num2str(tiempoPDZ3)])
