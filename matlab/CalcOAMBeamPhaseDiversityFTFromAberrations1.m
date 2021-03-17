%%
% File:	CalcOAMBeamPhaseDiversityFTFromAberrations1.m
% Author: Santiago Echeverri, Ren� Restrepo, Carlos Cuartas and N�stor
% Uribe
% Date:	  19/07/2016
% Modification:	  10/07/2020
% Vortex-enhanced coherent-illumination phase diversity for phase retrieval in coherent imaging systems 
% Optics Letters Vol. 41, Issue 8, pp. 1817-1820 (2016). https://doi.org/10.1364/OL.41.001817
%%
%  This program creates handles for every image in P.D. algorithm
%%
%  IN
%	FTSize: Size after creating a padarray with before fft's are computed.
%	beamDiameter: Pixels used to create grids.
%	gaussianC: Gaussian distribution waist.
%	rLimit: Maximum radius in Fourier field.
%	baseZAberrations: Vector array with base aberrations for PDs.
%	baseOAM: OAM guess (if unknow, otherwise = 0).
%	dx: Vortex displacement in x direction.
%	dy: Vortex displacement in y direction.
%	ZPD: Array with zernikes in each P.D.
%	OAMPD:  Vector with vortices topological charges.
%	nPDs: Number of Phase Diversities.
%
%  OUT
%	amplitudeCenter: Final amplitude.
%	IntensityCenter: Final intensity (|A^2|).
%	ampIni: Initial amplitude at object plane.
%	waveAberation: Surface with final aberrations.
%	gridSize: Simulated images size.
%	colRangeSLM: Column range at SLM.
%	rowRangeSLM: Row range at SLM.
%
% octave-3.8.1.
%%

%% Function definition

function [amplitudeCenter, IntensityCenter, ampIni, waveAberation, gridSize, colRangeSLM, rowRangeSLM] = CalcOAMBeamPhaseDiversityFTFromAberrations1(FTSize, beamDiameter, gaussianC, rLimit, baseZAberrations, baseOAM, dx, dy, ZPD, OAMPD, nPDs, coeff)

%% Check inputs

 size_baseZAberrations_2 = size(baseZAberrations, 2) ;
 size_baseZAberrations_1 = size(baseZAberrations, 1) ;

 %Check aberrations
  if size_baseZAberrations_2 > size_baseZAberrations_1
    baseZAberrations = baseZAberrations';
  end

 %Check Zernike PDs
  for k = 1:nPDs
    if size(ZPD{k}, 2) > size(ZPD{k}, 1)
      ZPD{k} = ZPD{k}';
    end
  end

 %Check OAM in PDs
  if size(OAMPD, 2) > size(OAMPD, 1)
    OAMPD = OAMPD';
  end

%% Create base images for each PD
 
  for k = 1:nPDs
    ZAberrations(:, k) = baseZAberrations + padarray(ZPD{k}, size_baseZAberrations_1- size(ZPD{k}, 1), 0, 'Post'); % Extend vector with initial aberrations so that they have same size as vector to retrieve

    [amplitudeCenter(:, :, k), IntensityCenter(:, :, k), ampIni, waveAberation(:, :, k), gridSize, colRangeSLM, rowRangeSLM] = CalcOAMBeamFTFromAberrations2(FTSize, beamDiameter, gaussianC, rLimit, ZAberrations(:, k), dx, dy, baseOAM + OAMPD(k), coeff); % Calculate images

    % --------------------------------------------

  end


  %% End

