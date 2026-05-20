Im = double(imread('pout.tif'));
    % Gradient image
h = fspecial('prewitt');
GImX = imfilter(Im,h,'symmetric');
GImY = imfilter(Im,h','symmetric');
GIm = sqrt(GImX.^2 + GImY.^2);
GIm = GIm/max(GIm(:));
    % Histogram
[Counts,Xh]=imhist(GIm);
    % Find max
cMax = find(Counts==max(Counts),1,'first');
Counts2 = Counts(cMax:end);
Xh2 = Xh(cMax:end);
    % Find threshold
tic;TPoint = TPointMethod(Counts2,Xh2,1),toc
tic;TPoint = TPointThreshold(Counts2,Xh2,2),toc
figure; plot(Xh,Counts), hold on, plot([Xh2(TPoint) Xh2(TPoint)],[0 max(Counts)],'r')
figure; 
    subplot(1,3,1), imshow(Im,[]), title('Initiale image')
    subplot(1,3,2), imshow(GIm,[]), title('Gradient image')
    subplot(1,3,3), imshow(GIm>Xh2(TPoint),[]), title('Segmented gradient image')
