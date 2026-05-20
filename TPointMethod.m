%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Created:  05/2008   Nicolas Coudray                                       %
%         Last modification:   1/2010  Nicolas Coudray                              %
%         Version 1.0                                                               %
%                                                                                   %
%         Copyright: 2006                                                           %
%         AnImATED TEM toolbox for HT3DEM   www.ht3dem.org                          %
%         Université de Haute Alsace - MIPS Laboratory - www.trop.mips.uha.fr       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function TPoint=TPointMethod(Counts,Xh,method)
% Searches the T-point of the descending or ascending curve Counts=f(Xh) with the 2 regression lines method:
% reference:
    % N. Coudray, J.-L. Buessler, J.-P. Urban
    % "Robust threshold estimation for images with unimodal histograms"
    % to be publised in "Pattern Recognition Letters" (2010)

% Inputs: 
%  - Counts, Xh: ordinates and abscissa of the ascending or descending part of the histogram;
%  - Method: if 1: use Matlab code
%               2: use faster compiled C code
% Output:
%  - TPoint: indice (in Xh) of the flection point 

if nargin==2
    method = 1;
elseif nargin~=3
    TPoint = -1;
    display('Error: invalid number of inputs')
    return
end
if length(Xh)~=length(Counts)
    TPoint = -1;
    display('Error: ordinates (Counts) and abscissa (Xh) must have the same length!')
    return
end

if method == 1
    nInd=find(Counts==max(Counts));
    if nInd>2   % the histogramm does not begin or end with a max
        if nInd<(length(Counts)-1)
            display('!! check max !!')
        end
    end

    if length(Xh)<=2 %check if the size of the data input is long enough
        TPoint=find(Xh==min(Xh),1,'first');
        return
    end

    for ii=2:length(Xh)-1       % for each abscissa...
        p1=SimplifyPolyfit(Xh(1:ii),Counts(1:ii),1);    %...calculate a linear regression with the points before....
        p2=SimplifyPolyfit(Xh(ii+1:end),Counts(ii+1:end),1);%... and another with the point after the abscisa...
        y1 = polyval(p1,Xh(1:ii));
        y2 = polyval(p2,Xh(ii+1:end));
        Erreur1(ii-1)=sum((y1-Counts(1:ii)).^2);                   %... mesure the square error for the two fittings...
        Erreur2(ii-1)=sum((y2-Counts(ii+1:end)).^2);
        ErreurFit(ii-1)=Erreur1(ii-1)+Erreur2(ii-1);      %... then add the errors
    end



    % Sum the errors
    TPoint=find(ErreurFit==min(ErreurFit),1,'first')+1;        % The T-Point is when the error for the two regression lines is minimal
else
    TPoint=TPointThreshold(Counts,Xh);
    if length(TPoint)>1
        TPoint=TPoint(1);
    end
end







function p=SimplifyPolyfit(x,y,n)   % simplification of Matlab polyfit
  
x = x(:);
y = y(:);
% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end
% Solve least squares problem.
[Q,R] = qr(V,0);
%ws = warning('off','all'); 
p = R\(Q'*y);    % Same as p = V\y;
%warning(ws);
r = y - V*p;
p = p.';          % Polynomial coefficients are row vectors by convention.
