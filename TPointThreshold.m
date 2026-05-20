%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Created:  2008   Jean-Luc Buessler                                        %
%         Last modification:  01/2008                                               %
%         Version 0.0                                                               %
%                                                                                   %
%         Copyright: 2006                                                           %
%         AnImATED TEM toolbox for HT3DEM   www.ht3dem.org                          %
%         Université de Haute Alsace - MIPS Laboratory - www.trop.mips.uha.fr       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TPointThreshold  determines the threshold in unimodal histogram using T_Point
%   
%   s = TPointThreshold(H)  where H is vector of bins from descending part
%   of an histogram equally spaced, returns indice of threshold
%
%   s = TPointThreshold(H, X) X is the vector of center values of the bins
%
%   [s ECoast]=TPointThreshold(H,...) also returns evaluation of the coast
%   function evaluated for each bin k (thresold is the arg min of ECoast)

%   Copyright 2008 University of Haute Alsace - HT3DEM project

% reference:
    % N. Coudray, J.-L. Buessler, J.-P. Urban
    % "Robust threshold estimation for images with unimodal histograms"
    % to be publised in "Pattern Recognition Letters" (2010)
