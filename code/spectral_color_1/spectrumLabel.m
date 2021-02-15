function varargout = spectrumLabel(hAxesTarget)
%spectrumLabel   Add a spectrum to the X-Axis of a plot.
%
%    H = SPECTRUMLABEL(HTARGET) adds a spectrum colorbar and wavelength
%    labels to the X axis of the axes object HTARGET.  For plots whose X
%    values are in the range 300 to 830, the colorbar will show the ROYGBIV
%    spectrum based on the 1964 CIE 10-degree observer.  Values outside the
%    range will be shown as black.  H contains a handle to the axes object
%    containing the spectrum.
%
%    Examples:
%    ---------
%    % Show the D65 illuminant.
%    [lambda, D65] = illuminant('d65');
%    figure;
%    plot(lambda, D65)
%    title('D_{65} illuminant')
%    spectrumLabel(gca)
%
%    % Show the CIE 1964 10-degree observer color matching functions.
%    [lambda, x, y, z] = colorMatchFcn('1964_full');
%    figure;
%    plot(lambda, x, 'r', lambda, y, 'g', lambda, z, 'b')
%    title('CIE 1964 10-degree observer')
%    spectrumLabel(gca)
%
%    See also colorMatchFcn, illuminant.

%    Copyright 1993-2005 The MathWorks, Inc.


hFig = ancestor(hAxesTarget, 'figure');

% Add a new axes to the same figure as the target axes.
figure(hFig);
hAxesSpectrum = axes;
set(hAxesSpectrum, 'visible', 'off')

% Areas outside the visible spectrum are black.
set(hAxesSpectrum, 'color', [1 1 1])

% Position the axes as appropriate.
targetPosition = get(hAxesTarget, 'position');
spectrumPosition = [targetPosition(1), ...
                    targetPosition(2), targetPosition(3), 1];
set(hAxesSpectrum, 'position', spectrumPosition)
set(hAxesSpectrum, 'units', 'pixels')

spectrumPosition = get(hAxesSpectrum, 'position');
% set(hAxesSpectrum, 'position', [spectrumPosition(1), ...                      % colorbar on bottom
%                                 spectrumPosition(2) - 20, ...
%                                 spectrumPosition(3), ...
%                                 20])

set(hAxesSpectrum, 'position', [spectrumPosition(1)+spectrumPosition(3)+25, ...    %color bar on right side
                                spectrumPosition(2)+30, ...
                                20, ...
                                spectrumPosition(4)-spectrumPosition(2)-90])
                            
% Line the X limits of the two axes up and display the spectrum.
set(hAxesSpectrum, 'ytick', get(hAxesTarget, 'ytick'))

[lambda, RGB] = createSpectrum('1964_full');
lambda = lambda(lambda>=400 & lambda<=700);     %show only visible range
RGB = RGB(:,41:341,:);                         %show only visible range
axes(hAxesSpectrum);
Eg = 1240./lambda;
frac = (length(Eg):-1:1)./length(Eg);
RGB1 = permute(RGB(1,:,:),[3 2 1]);
compare = [Eg; frac; RGB1]';
image(1:size(RGB,1), Eg, rot90(RGB));

Eg1 = fliplr(Eg(1:50:end));
Eg1 = round(Eg1,2);

% Turn off the unneeded axes labels.
set(hAxesSpectrum, 'xtick', [],'TickLength',[0 0]);
% set(hAxesSpectrum, 'yticklabels', flipud(get(hAxesSpectrum,'yticklabels'))) % color bar values on
% set(hAxesSpectrum, 'yticklabels', fliplr(get(gca,'ytick'))) % color bar values linear
% set(hAxesSpectrum, 'yticklabels',[]) % color bar values off
ticks = 3.1 - [0.04 0.2691 0.456 0.613 0.745 0.858 0.955].*(3.1-1.77);
set(hAxesSpectrum,'YTick',fliplr(ticks))
OppTickLabels = {'1.8' '2.0' '2.2' '2.4' '2.6' '2.8' '3.0'};
set(hAxesSpectrum,'YTickLabel',fliplr(OppTickLabels),'FontSize',15);
set(hAxesSpectrum,'YAxisLocation','right')

% Make the figure visible.
set(hAxesSpectrum, 'units', 'normalized')
set(hAxesSpectrum, 'visible', 'on')

% Return the handle if requested.
if (nargout == 1)
    
    varargout{1} = hAxesSpectrum;
    
end

% Add a callback to update the label if the X limits of the target changes.
% TO DO.