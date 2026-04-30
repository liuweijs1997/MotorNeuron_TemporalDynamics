function gdspk = selectspike(spk)
% SELECTSPIKE - Interactive posthoc spike classification
%    gdspk = SELECTSPIKE(spk) plots a raster of the spikes in SPK and lets
%    the user interactively select which spikes belong to the neuron of
%    interest.
%
%    Place the green and blue lines around the relevant dots (it does not
%    matter which one is above and which one is below). More handles
%    can be made by dragging the line; handles can be removed by dragging
%    them past the next handle.
%
% Updated for newer MATLAB versions:
% - do not use figure handle as cell index
% - replace old property abbreviations with full property names

global phsc_data

if isempty(phsc_data)
    phsc_data = {};
end

[w,h] = screensize;
figw = w - 40;
figh = h - 80;

f = figure( ...
    'Units', 'pixels', ...
    'Position', [20 40 figw figh], ...
    'MenuBar', 'none', ...
    'NumberTitle', 'off', ...
    'Name', 'Posthoc_SpikeClass');

% use numeric index for phsc_data
idx = numel(phsc_data) + 1;
setappdata(f, 'phsc_idx', idx);

phsc_data{idx}.loaded = 0;
phsc_data{idx}.figw = figw;
phsc_data{idx}.figh = figh;
phsc_data{idx}.figure = f;

butw = 90;

uicontrol( ...
    'String', 'Done', ...
    'Position', [10 figh-30 butw 25], ...
    'Style', 'pushbutton', ...
    'Callback', @phsc_done);

uicontrol( ...
    'String', 'Zoom', ...
    'Position', [figw-butw-10 figh-30 butw 25], ...
    'Tag', 'zoom', ...
    'Style', 'checkbox', ...
    'Value', 1);

phsc_data{idx}.axesw = figw - 60;
phsc_data{idx}.axesh = figh - 90;

axes( ...
    'Units', 'pixels', ...
    'Position', [55 45 phsc_data{idx}.axesw phsc_data{idx}.axesh], ...
    'Tag', 'graph', ...
    'ButtonDownFcn', @phsc_click);

hold on
imagesc(zeros(10,10), 'Tag', 'trace', 'ButtonDownFcn', @phsc_click);
plot([0 0],[0 0], 'b', 'LineWidth', 2, 'Tag', 'lowerline', 'ButtonDownFcn', @phsc_click);
plot(0,0, 'b.', 'MarkerSize', 20, 'Tag', 'lowerdots', 'ButtonDownFcn', @phsc_click);
plot([0 0],[0 0], 'g', 'LineWidth', 2, 'Tag', 'upperline', 'ButtonDownFcn', @phsc_click);
plot(0,0, 'g.', 'MarkerSize', 20, 'Tag', 'upperdots', 'ButtonDownFcn', @phsc_click);

a = axis;
setappdata(f, 'axlim0', a);
setappdata(f, 'axlim', a);

phsc_loaddat(f, spk);
setappdata(f, 'completed', 0);

while true
    uiwait(f);
    cancel = 0;
    done = 0;
    try
        done = getappdata(f, 'completed');
    catch
        cancel = 1;
    end
    if done || cancel
        break
    end
end

if done
    gdspk = phsc_getdata(f);
    close(f);
else
    gdspk.tms = [];
    gdspk.amp = [];
end

end


%----------------------------------------------------------------------
function phsc_loaddat(figh, spk)
global phsc_data
idx = getappdata(figh, 'phsc_idx');

phsc_data{idx}.ifn = 'data';
phsc_data{idx}.src.spk.tms = spk.tms;
phsc_data{idx}.src.spk.chs = 1 + 0 * spk.tms;
phsc_data{idx}.src.spk.hei = spk.amp;
phsc_data{idx}.src.type = 'spikes';

phsc_process_spikes(figh);
end


%----------------------------------------------------------------------
function phsc_process_spikes(figh)
global phsc_data
idx = getappdata(figh, 'phsc_idx');

phsc_data{idx}.loaded = 1;
phsc_data{idx}.c = 1;

switch phsc_data{idx}.src.type
    case 'spikes'
        phsc_data{idx}.C = 1;
        phsc_data{idx}.T = max([max(phsc_data{idx}.src.spk.tms) 1]);
    case 'histo'
        phsc_data{idx}.C = size(phsc_data{idx}.src.hst, 3);
        phsc_data{idx}.T = size(phsc_data{idx}.src.hst, 2);
end

phsc_data{idx}.lower_thr = cell(1, phsc_data{idx}.C);
phsc_data{idx}.upper_thr = cell(1, phsc_data{idx}.C);

if isfield(phsc_data{idx}.src, 'chnames')
    phsc_data{idx}.chnames = phsc_data{idx}.src.chnames;
else
    phsc_data{idx}.chnames = cell(1, phsc_data{idx}.C);
end


dT = min(60, phsc_data{idx}.T / 2);

yyAll = phsc_data{idx}.src.spk.hei(:);
yyAll = yyAll(~isnan(yyAll));

for c = 1:phsc_data{idx}.C
    % robust initialization of threshold lines
    low0 = prctile(yyAll, 5);
    high0 = prctile(yyAll, 95);

    phsc_data{idx}.lower_thr{c} = [dT low0];
    phsc_data{idx}.upper_thr{c} = [dT high0];
end

delete(findobj(figh, 'Tag', 'channelselect'));

figure(figh);
if phsc_data{idx}.C > 1
    for c = 1:phsc_data{idx}.C
        h = uicontrol( ...
            'String', sprintf('Ch%i %s', c, phsc_data{idx}.chnames{c}), ...
            'Style', 'pushbutton', ...
            'Tag', 'channelselect', ...
            'UserData', c, ...
            'Position', [500 + c*75 phsc_data{idx}.figh-30 70 25], ...
            'Callback', @phsc_channelselect);
        if c == phsc_data{idx}.c
            set(h, 'FontWeight', 'bold');
        end
    end
end

set(figh, 'Name', sprintf('Posthoc_SpikeClass: %s', phsc_data{idx}.ifn));

minY = min(yyAll) - 5;
maxY = max(yyAll) + 5;

setappdata(figh, 'axlim0', [0 phsc_data{idx}.T minY maxY]);
setappdata(figh, 'axlim',  [0 phsc_data{idx}.T minY maxY]);

phsc_redraw(figh, 1);
end


%----------------------------------------------------------------------
function phsc_redraw(figh, graphtoo)
global phsc_data
idx = getappdata(figh, 'phsc_idx');

axes(findobj(figh, 'Tag', 'graph'));
a = getappdata(figh, 'axlim');
axis(a);

if graphtoo
    c = phsc_data{idx}.c;

    switch phsc_data{idx}.src.type
        case 'spikes'
            ii = find(phsc_data{idx}.src.spk.chs == c);
            xx = phsc_data{idx}.src.spk.tms(ii);
            yy = phsc_data{idx}.src.spk.hei(ii);

            a = axis;
            ii = find(xx >= a(1) & xx <= a(2) & yy >= a(3) & yy <= a(4));

            if isempty(ii)
                nn = zeros(10,10);
                xx = linspace(a(1), a(2), 10);
                yy = linspace(a(3), a(4), 10);
            else
                [nn, xx, yy] = hist2(xx(ii), yy(ii), phsc_data{idx}.axesw, phsc_data{idx}.axesh);
                nn = gsmooth(gsmooth(nn, 2.5)', 1)';
            end

        case 'histo'
            yy = sw_sig2log(phsc_data{idx}.src.hst_xx);
            xx = 1:phsc_data{idx}.T;
            nn = phsc_data{idx}.src.hst(:, :, c);

            a = axis;
            xidx = find(xx >= a(1) & xx <= a(2));
            yidx = find(yy >= a(3) & yy <= a(4));
            xx = xx(xidx);
            yy = yy(yidx);
            nn = nn(yidx, xidx);

            scl = ceil(length(xx) / phsc_data{idx}.axesw);
            bin = floor(length(xx) / scl);

            if scl * bin > 0
                xx = mean(reshape(xx(1:scl*bin), [scl bin]), 1);
                nn = squeeze(sum(reshape(nn(:,1:scl*bin), [length(yy) scl bin]), 2));
                nn = gsmooth(gsmooth(nn, .05)', .5)';
            else
                nn = zeros(10,10);
                xx = linspace(a(1), a(2), 10);
                yy = linspace(a(3), a(4), 10);
            end
    end

    h = findobj(figh, 'Tag', 'trace');
    set(h, 'XData', xx, 'YData', yy, 'CData', nn);

    xlabel('Time (sec)');
    ylabel('Spike Ampl.');
    colormap(hotpow(200, .25));
    
    nn1 = sort(nn(:));
    if ~isempty(nn1)
        caxis([0 max(1, nn1(ceil(length(nn1) * .975)))]);
    end
    
    set(gca, ...
        'TickDir', 'out', ...
        'TickLength', [.004 .002]);
    ytickformat('%.0f');
    
    axis(a);
end

a = axis;
c = phsc_data{idx}.c;

% lower line
xx = phsc_data{idx}.lower_thr{c}(:,1);
yy = phsc_data{idx}.lower_thr{c}(:,2);

h = findobj(gca, 'Tag', 'lowerdots');
set(h, 'XData', xx, 'YData', yy);

h = findobj(gca, 'Tag', 'lowerline');
set(h, 'XData', [a(1); xx; a(2)], 'YData', [yy(1); yy; yy(end)]);

% upper line
xx = phsc_data{idx}.upper_thr{c}(:,1);
yy = phsc_data{idx}.upper_thr{c}(:,2);

h = findobj(gca, 'Tag', 'upperdots');
set(h, 'XData', xx, 'YData', yy);

h = findobj(gca, 'Tag', 'upperline');
set(h, 'XData', [a(1); xx; a(2)], 'YData', [yy(1); yy; yy(end)]);
end


%----------------------------------------------------------------------
function phsc_click(h, x) %#ok<INUSD>
global phsc_data
figh = gcbf;
idx = getappdata(figh, 'phsc_idx');

iszoom = get(findobj(figh, 'Tag', 'zoom'), 'Value');
c = phsc_data{idx}.c;
tag = get(h, 'Tag');
act = 0;

switch tag
    case 'lowerdots'
        xx = phsc_data{idx}.lower_thr{c}(:,1);
        yy = phsc_data{idx}.lower_thr{c}(:,2);
        xy = get(gca, 'CurrentPoint'); 
        xy = xy(1,1:2);
        [~, ii] = min((xy(1)-xx).^2 + (xy(2)-yy).^2);
        mousemove(@phsc_move, 'lower', xx, yy, ii);
        act = 1;

    case 'lowerline'
        xx = phsc_data{idx}.lower_thr{c}(:,1);
        yy = phsc_data{idx}.lower_thr{c}(:,2);
        xy = get(gca, 'CurrentPoint'); 
        xy = xy(1,1:2);
        ii = findfirst_ge(xx, xy(1));
        if ii > 0
            xx = [xx(1:ii-1); xy(1); xx(ii:end)];
            yy = [yy(1:ii-1); xy(2); yy(ii:end)];
        else
            xx = [xx; xy(1)];
            yy = [yy; xy(2)];
            ii = length(xx);
        end
        mousemove(@phsc_move, 'lower', xx, yy, ii);
        act = 1;

    case 'upperdots'
        xx = phsc_data{idx}.upper_thr{c}(:,1);
        yy = phsc_data{idx}.upper_thr{c}(:,2);
        xy = get(gca, 'CurrentPoint'); 
        xy = xy(1,1:2);
        [~, ii] = min((xy(1)-xx).^2 + (xy(2)-yy).^2);
        mousemove(@phsc_move, 'upper', xx, yy, ii);
        act = 1;

    case 'upperline'
        xx = phsc_data{idx}.upper_thr{c}(:,1);
        yy = phsc_data{idx}.upper_thr{c}(:,2);
        xy = get(gca, 'CurrentPoint'); 
        xy = xy(1,1:2);
        ii = findfirst_ge(xx, xy(1));
        if ii > 0
            xx = [xx(1:ii-1); xy(1); xx(ii:end)];
            yy = [yy(1:ii-1); xy(2); yy(ii:end)];
        else
            xx = [xx; xy(1)];
            yy = [yy; xy(2)];
            ii = length(xx);
        end
        mousemove(@phsc_move, 'upper', xx, yy, ii);
        act = 1;

    otherwise
        if iszoom
            switch get(gcbf, 'SelectionType')
                case {'alt', 'open'}
                    setappdata(figh, 'axlim', getappdata(figh, 'axlim0'));
                otherwise
                    xy = get(gca, 'CurrentPoint'); 
                    xy = xy(1,1:2);
                    rbbox;
                    xy1 = get(gca, 'CurrentPoint'); 
                    xy1 = xy1(1,1:2);
                    dd = (xy1(1)-xy(1))^2 + (xy1(2)-xy(2))^2;
                    if dd > 5
                        setappdata(figh, 'axlim', [sort([xy(1) xy1(1)]) sort([xy(2) xy1(2)])]);
                    else
                        setappdata(figh, 'axlim', getappdata(figh, 'axlim0'));
                    end
            end
        end
        phsc_redraw(figh, 1);
end

if act
    xx = phsc_data{idx}.upper_thr{c}(:,1);
    yy = phsc_data{idx}.upper_thr{c}(:,2);
    nn = find(xx(2:end) < xx(1:end-1));
    xx(nn) = (xx(nn+1)+xx(nn))/2; 
    xx(nn+1) = [];
    yy(nn) = (yy(nn+1)+yy(nn))/2; 
    yy(nn+1) = [];
    phsc_data{idx}.upper_thr{c} = [xx yy];

    xx = phsc_data{idx}.lower_thr{c}(:,1);
    yy = phsc_data{idx}.lower_thr{c}(:,2);
    nn = find(xx(2:end) < xx(1:end-1));
    xx(nn) = (xx(nn+1)+xx(nn))/2; 
    xx(nn+1) = [];
    yy(nn) = (yy(nn+1)+yy(nn))/2; 
    yy(nn+1) = [];
    phsc_data{idx}.lower_thr{c} = [xx yy];

    phsc_redraw(figh, 0);
end
end


%----------------------------------------------------------------------
function phsc_move(h, xy0, xy1, whch, xx, yy, ii) %#ok<INUSD>
global phsc_data
figh = get(h, 'Parent');
idx = getappdata(figh, 'phsc_idx');
c = phsc_data{idx}.c;

xx(ii) = xy1(1);
yy(ii) = xy1(2);

switch whch
    case 'lower'
        phsc_data{idx}.lower_thr{c} = [xx yy];
    case 'upper'
        phsc_data{idx}.upper_thr{c} = [xx yy];
end

phsc_redraw(figh, 0);
end


%----------------------------------------------------------------------
function phsc_channelselect(h, x) %#ok<INUSD>
global phsc_data
figh = gcbf;
idx = getappdata(figh, 'phsc_idx');

phsc_data{idx}.c = get(gcbo, 'UserData');
set(findobj(gcbf, 'Tag', 'channelselect'), 'FontWeight', 'normal');
set(gcbo, 'FontWeight', 'bold');
phsc_redraw(figh, 1);
end


%----------------------------------------------------------------------
function phsc_done(h, x) %#ok<INUSD>
setappdata(gcbf, 'completed', 1);
uiresume(gcbf);
end


%----------------------------------------------------------------------
function gdspk = phsc_getdata(figh)
global phsc_data
idx = getappdata(figh, 'phsc_idx');

gdspk.tms = [];
gdspk.amp = [];

if strcmp(phsc_data{idx}.src.type, 'spikes')
    for c = 1:phsc_data{idx}.C
        ii = find(phsc_data{idx}.src.spk.chs == c);
        tt = phsc_data{idx}.src.spk.tms(ii);
        yy = phsc_data{idx}.src.spk.hei(ii);

        x = phsc_data{idx}.lower_thr{c}(:,1);
        y = phsc_data{idx}.lower_thr{c}(:,2);
        [x, ord] = sort(x);
        y = y(ord);
        x = [min([0 min(x)-10]); x(:); max([max(tt) max(x)]) + 10];
        y = [y(1); y(:); y(end)];
        lower_thr = interp1(x, y, tt, 'linear');

        x = phsc_data{idx}.upper_thr{c}(:,1);
        y = phsc_data{idx}.upper_thr{c}(:,2);
        [x, ord] = sort(x);
        y = y(ord);
        x = [min([0 min(x)-1]); x(:); max([max(tt) max(x)]) + 10];
        y = [y(1); y(:); y(end)];
        upper_thr = interp1(x, y, tt, 'linear');

        if mean(upper_thr) < mean(lower_thr)
            ll = lower_thr;
            lower_thr = upper_thr;
            upper_thr = ll;
        end

        ii = ii(yy > lower_thr & yy < upper_thr);

        gdspk.tms = phsc_data{idx}.src.spk.tms(ii);
        gdspk.amp = phsc_data{idx}.src.spk.hei(ii);
    end
end
end