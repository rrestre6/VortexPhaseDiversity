%%
% File:	  FTBeamPhaseDiversityErrorFunction1.m
% Author: Santiago Echeverri, René Restrepo, Carlos Cuartas and Néstor
% Uribe
% Date:	  19/07/2016
% Modification:	  10/07/2020
% Vortex-enhanced coherent-illumination phase diversity for phase retrieval in coherent imaging systems 
% Optics Letters Vol. 41, Issue 8, pp. 1817-1820 (2016). https://doi.org/10.1364/OL.41.001817
%%
%  This program calculates image difference (error) in each iteration
%%
%  IN
%	FTSize: Size after creating a padarray with before fft's are computed.
%	beamDiameter: Pixels used to create grids.
%	gaussianC: Gaussian distribution waist.
%	rLimit: Maximum radius in Fourier field.
%	baseZAberrationsDxDy: Vector array with base aberrations for PDs.
%	baseOAM: OAM guess (if unknow, otherwise = 0).
%	ZPD: Array with zernikes in each P.D.
%	OAMPD:  Vector with vortices topological charges.
%	nPDs: Number of Phase Diversities.
%	ExpIntensity: Experimental images array.
%	colRangeIn: Column range at the entrace field.
%	rowRangeIn: Row range at the entrace field.
%	hAbsFT: If verbosity, images to change in each iteration showing intensity.
%	hPhaseIFT: If verbosity, image to change in each iteration showing phase.
%	fitDxDy: if known dx and dy this is not calculated.
%	verbosity: Verbosity level.
%
%  OUT
%	error: Numerical error between simulated an real intensities.
%
% octave-3.8.1.
%%

%% Function definition

function [error] = FTBeamPhaseDiversityErrorFunction1(FTSize, beamDiameter, gaussianC, rLimit, baseZAberrationsDxDy, baseOAM, ZPD, OAMPD, nPDs, ExpIntensity, colRangeIn, rowRangeIn, hAbsFT, hPhaseIFT, fitDxDy, verbosity)

%% Check if dx and dy are already known

  if fitDxDy
    baseZAberrations = baseZAberrationsDxDy(1:(end - 2));
    dx = baseZAberrationsDxDy(end - 1);
    dy = baseZAberrationsDxDy(end);
  else
    baseZAberrations = baseZAberrationsDxDy;
    dx = 0; dy = 0;
  end

  %% Create simulated images after changing baseZAberrations to fit the best guest
  [amplitudeCenter, IntensityCenter, ampIni, waveAberation, gridSize, colRangeSLM, rowRangeSLM] = CalcOAMBeamPhaseDiversityFTFromAberrations1(FTSize, beamDiameter, gaussianC, rLimit, baseZAberrations, baseOAM, dx, dy, ZPD, OAMPD, nPDs);

  %% Error calculus
%The error is taken from the square sum between all the simulated
%and all the experimental images. The minimization search the
%minimun error. 
  
  nmPixels =  size(ExpIntensity(colRangeIn, rowRangeIn, 1) ,1) * ...
      size(ExpIntensity(colRangeIn, rowRangeIn, 1) ,2);
  
  error = 0;
  for k = 1:nPDs
    error = error + sum(sum(abs(ExpIntensity(colRangeIn, rowRangeIn, ...
                                              k) - ...
                                 IntensityCenter(colRangeIn, rowRangeIn, k)) .^ 2, 1), 2) /nmPixels;  
    if verbosity
      set(hAbsFT(k), 'cdata', abs(IntensityCenter(colRangeIn, rowRangeIn, k)));
    end
  end
  if verbosity
    set(hPhaseIFT, 'cdata', (angle(waveAberation(:,:,1))));
    drawnow;
  end

