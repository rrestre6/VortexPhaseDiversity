%%
% File: demo_RetrievePhaseFromPhaseDiversity1.m
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
%		the procedure is employed using the two best-fitting diversities $-2, 2$
%		The idea is to get 1lambda astigmatism
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
rLimit = 74/gridSize;% FT gaussian effective size
verbosity = 1;
seed = 0;
dx = 0;
dy = 0;
%mType = 'nelder';
mType = 'nonlin';
mPars = {1000 1e-12 4000 300};% Minimization parameters, for this procedure [isz, tolFun]
% The isz is the initial step size which determines the magnitude of the steps in first iteration - recommended 1-
% tolFun is the tolerance of the function determined by the change between two successive iteration, if the error change is
% smaller than the given tolerance, the algorithm stops. -recommended range [1e-4 - 1e-7]-.


%% Diversities
% As we only want to use topological charge as diversity, $l$ diversities are [-2, 2], then $k$ diversities are vectors inside a cell
% which contains the diversities vectors, zero in this case

% For lDiv = [-2, 0, 2] uncomment. 
kDiv = {[0 0 0 0 0 0 0 0 0 0 0 0 0]', [0 0 0 0 0 0 0 0 0 0 0 0 0]', [0 0 0 0 0 0 0 0 0 0 0 0 0]'};
lDiv = [-2, 0, 2];

% For lDiv = [-2, -1, 0, 1, 2] uncomment
% kDiv = {[0 0 0 0 0 0 0 0 0 0 0 0 0]', [0 0 0 0 0 0 0 0 0 0 0 0 0]', [0 0 0 0 0 0 0 0 0 0 0 0 0]', ...
%    [0 0 0 0 0 0 0 0 0 0 0 0 0]',[0 0 0 0 0 0 0 0 0 0 0 0 0]'};
% lDiv = [-2, -1, 0, 1, 2];

%% Desired aberration
% In our notation, aberrations follow
% [defocus, astigmatism, astigmatism45°, comaX, comaY, spherical, trefoil, trefoil45°, ...
% secondAstigmatism, secondAstigmatism45°, secondComaX, secondComaY, ...],
ZAberration = [0 0 0 0 0 0 0 0 0 0 0 0 0]'; % This is astigmatism

%% Image range
colRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;
rowRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;

%% Diversity images creation
% We use CalcOAMBeamFTFromAberrations2 which is the image simulator in order to produce the desire diversities.
% This one has OAM -2

addpath(genpath('../CoherentPhaseDiversity_DistributionMatlabVersion'));
 
% TC0 = double(fread(fopen('./OAMs/20210303/TC0.bin'), [256 256], 'uint8'));
% TC1 = double(fread(fopen('./OAMs/20210303/TC1.bin'), [256 256], 'uint8'));
% TC_1 = double(fread(fopen('./OAMs/20210303/TC-1.bin'), [256 256], 'uint8'));
% TC2 = double(fread(fopen('./OAMs/20210303/TC2.bin'), [256 256], 'uint8'));
% TC_2 = double(fread(fopen('./OAMs/20210303/TC-2.bin'), [256 256], 'uint8'));

% TC0 = double(fread(fopen('./OAMs/20210311_TwoOptimizations/TC0_CC.bin'), [256 256], 'uint8'));
% TC1 = double(fread(fopen('./OAMs/20210311_TwoOptimizations/TC1_CC.bin'), [256 256], 'uint8'));
% TC_1 = double(fread(fopen('./OAMs/20210311_TwoOptimizations/TC-1_CC.bin'), [256 256], 'uint8'));
% TC2 = double(fread(fopen('./OAMs/20210311_TwoOptimizations/TC2_CC.bin'), [256 256], 'uint8'));
% TC_2 = double(fread(fopen('./OAMs/20210311_TwoOptimizations/TC-2_CC.bin'), [256 256], 'uint8'));

TC0 = double(fread(fopen('./OAMs/20210304_simulados/TC0.bin'), [120 120], 'double'));
TC1 = double(fread(fopen('./OAMs/20210304_simulados/TC1.bin'), [120 120], 'double'));
TC_1 = double(fread(fopen('./OAMs/20210304_simulados/TC-1.bin'), [120 120], 'double'));
TC2 = double(fread(fopen('./OAMs/20210304_simulados/TC2.bin'), [120 120], 'double'));
TC_2 = double(fread(fopen('./OAMs/20210304_simulados/TC-2.bin'), [120 120], 'double'));

[maxval,idx]=max(TC0(:));
[row,col]=ind2sub(size(TC0), idx);

xSize = [row-60:row+59];
ySize = [col-60:col+59];

% For lDiv = [-2, 0, 2] uncomment
ROIImageZ3 =zeros(120,120,3);
ROIImageZ3(:,:,1) = (TC_2(xSize,ySize))/(max(TC_2(:)));
ROIImageZ3(:,:,2) = (TC0(xSize,ySize))/(max(TC0(:)));
ROIImageZ3(:,:,3) = (TC2(xSize,ySize))/(max(TC2(:)));

% For lDiv = [-2, -1, 0, 1, 2] uncomment
% ROIImageZ3 =zeros(120,120,5);
% ROIImageZ3(:,:,1) = (TC_2(xSize,ySize))/(max(TC_2(:)));
% ROIImageZ3(:,:,2) = (TC_1(xSize,ySize))/(max(TC_1(:)));
% ROIImageZ3(:,:,3) = (TC0(xSize,ySize))/(max(TC0(:)));
% ROIImageZ3(:,:,4) = (TC1(xSize,ySize))/(max(TC1(:)));
% ROIImageZ3(:,:,5) = (TC2(xSize,ySize))/(max(TC2(:)));

%% PD run
% For experimental images replace ROIImageZ3 with your images organized as the diversities

% Initial seed for the algorithm
if seed == 1
    lDiversities = [-1 1];
    kDiversities = {[0 0 0]' [0 0 0 0]'};
    beamImage =zeros(120,120,2);
    beamImage(:,:,1) = (TC0(xSize,ySize))/(max(TC0(:)));
    beamImage(:,:,2) = (TC1(xSize,ySize))/(max(TC1(:)));
    [Phase, intensity, ZBaseAberration] = RetrievePhaseFromSpiralPDModified3(beamImage, FTSize,...
        gridSize, gaussianC, rLimit, lDiversities, kDiversities, mPars{4}, verbosity);
    ZBaseAberration = ZBaseAberration';
    
else
    ZBaseAberration = zeros(15,1); % Usually 15 zernikes zeros(15,1)
    ZBaseAberration(1)=5;  ZBaseAberration(7)=0.3;  ZBaseAberration(8)=0.4;
end

tic; % time
[ZAberrationsZ3, ~, ~, ~, ~, ~, ~, ~, waveAberationPD, PDIntensityZ3] = RetrievePhaseFromPhaseDiversity1(ROIImageZ3,...
    FTSize, gridSize, gaussianC, rLimit, ZBaseAberration, 0, dx, dy, ...
    kDiv, lDiv, mType, mPars, verbosity); % Run phase diversities

% note that printed in terminal is the final error and the aberrations retrieved, which clearly is 1lambda astigmatism
tiempoPDZ3 = toc;

% Final result info
errPDZ3 = sum(sum(abs(ROIImageZ3(:,:,1) - PDIntensityZ3(colRangeSLM,rowRangeSLM,2)) / (size(PDIntensityZ3,1) * size(PDIntensityZ3,2)), 1), 2);
disp(['Error final por PDZ3 ',num2str(errPDZ3),' tiempo ', num2str(tiempoPDZ3)])

