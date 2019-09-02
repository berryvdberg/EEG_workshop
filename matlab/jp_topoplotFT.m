function h = jp_topoplotFT(cfg,varargin)

% JP_TOPOPLOT_FT plots topographic maps using linear griddata
% interpolation. Based on Biharmonic spline interpolation as implemented in
% matlab
%
% Use as
% JP_TOPOPLOT(cfg,data,....)
%
% The input data should be organised in a structure as obtained from
% the FT_FREQGRANDAVERAGE or the FT_TIMELOCKGRANDAVERAGE function.
%
% The configuration should contain:
% cfg.latency = M*2 matrix with the time bins to be plotted: [0   0.1
%                                                             0.1 0.2
%                                                             0.1 0.2
%                                                             0.1 0.2]
%
% cfg.foi     = frequency of interest, only applicable if dimord contains
%               'freq'
% cfg.elecsize = size of the electrodes
% cfg.perspective = Mx1 cell-array with for each data input a perspective
%                   {'back' 'back'  'back'};
% cfg.ncontours = the amount of contour lines to be plotted
% cfg.granularity = how fine the grid should be, lower numbers equals nicer
%                   figures
% cfg.plotchannellab = 'yes' or 'no' whether or not to plot channel labels
% cfg.cutoff  = how low is the data plotted on the side views and back
%               view, a number between -1 and 1, or a channel label
% cfg.clim = the color scale limit
% cfg.layout = should contain the elec structure including a .pnt structure
%              with the xyz coordinates of each electrode
%
% Joonkoo Park (joonkoo@umass.edu)
% Berry van den Berg (berryv.dberg@gmail.com)

% order the data
if nargin > 1
    cfg.dataname = {inputname(2)};
    for k = 3:nargin
        cfg.dataname{end+1} = inputname(k);
    end
end
Ndata = numel(varargin);

% default granularity
if isfield(cfg,'granularity')
    meshgrn = cfg.granularity;    % too fine granularity will slow down the plotting
else
    meshgrn = 0.05;
end

grn = 0.005;

% number of topo maps
numtopo = size(cfg.latency,1);

% number of bins
numbin  = Ndata;
h = figure('Position',[100 100 200*numtopo 250*numbin]);

for ibin = 1 : numbin
    % get the xyz coordinates here
    xyz = zeros(length(varargin{ibin}.label),3);
    for i = 1:length(varargin{ibin}.label);
        xyz(i,:) = cfg.layout.pnt(strcmp(varargin{ibin}.label{i},cfg.layout.label),:);
    end

    % normalize
    chanlabel = varargin{ibin}.label;
    xyz = bsxfun(@rdivide,xyz,sqrt(sum(xyz.^2,2)));

    % start loop for the number of topomaps
    for itopo = 1 : numtopo
        subplot(numbin,numtopo,(ibin-1)*numtopo+itopo);

        % get time index for topo
        idxtime = varargin{ibin}.time >= cfg.latency(itopo,1) & varargin{ibin}.time <= cfg.latency(itopo,2);

        % prepare data 'v', here you can add mode data types for different
        % dimord
        if strcmp(varargin{ibin}.dimord, 'subj_chan_time');
            v = mean(mean(varargin{ibin}.individual(:,:,idxtime),1),3);
        elseif strcmp(varargin{ibin}.dimord, 'subj_chan_freq_time');
            idxfreq = varargin{ibin}.freq >= cfg.foi(1) & varargin{ibin}.freq <= cfg.foi(2);
            v = mean(mean(mean(varargin{ibin}.powspctrm(:,:,idxfreq,idxtime),1),3),4);
        elseif strcmp(varargin{ibin}.dimord, 'chan_freq_time');
            idxfreq = varargin{ibin}.freq >= cfg.foi(1) & varargin{ibin}.freq <= cfg.foi(2);
            v = mean(mean(varargin{ibin}.powspctrm(:,idxfreq,idxtime),2),3);
        elseif strcmp(varargin{ibin}.dimord, 'chan_time');
            v =mean(varargin{ibin}.individual(:,idxtime),2);
        else
            error('no support for dimord')
        end



        % start the topo plotting
        switch cfg.perspective{ibin}
            case 'back'

                % Back view is where X < 0
                chanidx = find(xyz(:,1)<=0.01);
                xyz_b = xyz(chanidx,:);
                chanlabelb = chanlabel(chanidx);
                x = -xyz_b(:,2); % *-1 is because you are looking from the back
                y = xyz_b(:,3);

                [xq,yq] = meshgrid(-1:meshgrn:1, -1:meshgrn:1);

                vq = griddata(x,y,v(chanidx),xq,yq,'v4');

                % trim outside the unit circle and the minimum y val
                vq((xq.^2 + yq.^2) > 1) = NaN;
                vq((.15 .* xq.^2 + yq + .6) < 0) = NaN;
                vq((.4 .* xq.^2 + yq + 0.3) < 0) = NaN;

                % define a cutoff where below the plotting should be white
      % define a cutoff where below the plotting should be white
                       % in case of 10-10 this should not be A1 or A2
                if isfield(cfg,'cutoff')
                    % if the cutoff is hardcoded
                    if isnumeric(cfg.cutoff)
                        cutoff = cfg.cutoff;
                    else
                        % if the cutoff is a channel label
                        cutoff = ismember(chanlabel(chanidx),cfg.cutoff);
                        cutoff = y(cutoff);
                    end
                    vq(yq < cutoff) = NaN;

                else
% 
%                     vq((.15 .* xq.^2 + yq + .6) < 0) = NaN;
%                              %       vq((xq.^2 + yq.^2) > 1) = NaN;

                end
                % plot the contour maps, without lines
                contourf(xq,yq,vq,cfg.ncontours,'fill','on','LineStyle','none');

%                 pcolor(xq(1,:),yq(:,1),vq);
%                 shading flat
%                 shading interp
% %
                if isempty(cfg.clim)
                    maxval = ceil(max(abs(v)) * 10) / 10;
                    caxis([-maxval maxval]);
                else
                    caxis(cfg.clim(:)); % add support for different clim per bin?
                end

                %colorbar;
                hold on;
                % scatter(-xyz_b(:,2),xyz_b(:,3),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');

                scatter(-xyz_b(:,2),xyz_b(:,3),cfg.elecsize,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');

                % head & neck (hn_)
                hn_k = 0.45;     % head-neck joint
                hn_b = 1.15;     % mid-point of neck
                hn_x1 = -1:grn:hn_k;
                hn_x2 = hn_k:grn:1.2;
                hn_y1 = sqrt(1-hn_x1.^2);
                hn_a = (-hn_k / sqrt(1-hn_k.^2)) / (2*(hn_k-hn_b));     % curve of parabola
                hn_c = sqrt(1-hn_k.^2) - hn_a .* (hn_k - hn_b) .^ 2;
                hn_y2 = hn_a * (hn_x2 - hn_b) .^ 2 + hn_c;


                hn_y = -[hn_x1, hn_x2];
                hn_x = [hn_y1, hn_y2];

%                 % ears
%                 ear_ytop = -0.2;
%                 ear_xtop = -hn_x(hn_y == ear_ytop);
%                 ear_ybottom = -0.8;
%                 ear_xbottom = -hn_x(hn_y == ear_ybottom);
%
%                 ear_y =  ear_ybottom:grn:ear_ytop;
%                 ear_x = linspace(ear_xbottom,ear_xtop,numel(ear_y));
%
%
%                 figure;
%                 r = 1;
%                 x=-r:.01:r;
%                 y= sqrt(r^2-x.^2);
%                 hold on;
%                 plot(x,y);
%                 % first element, bottom ear
%                 x1_start = -0.3
%
%                 y1 =
%                 x1 = -1:grn:hn_k;
%
%                                 % head & neck (hn_)
%                 hn_k = 0.45;     % head-neck joint
%                 hn_b = 1.15;     % mid-point of neck
%                 hn_x1 = -1:grn:hn_k;
%                 hn_x2 = hn_k:grn:1.2;
%                 hn_y1 = sqrt(1-hn_x1.^2);
%                 hn_a = (-hn_k / sqrt(1-hn_k.^2)) / (2*(hn_k-hn_b));     % curve of parabola
%                 hn_c = sqrt(1-hn_k.^2) - hn_a .* (hn_k - hn_b) .^ 2;
%                 hn_y2 = hn_a * (hn_x2 - hn_b) .^ 2 + hn_c;
%
%                 figure
%                 plot(hn_x1,hn_y1);


                %plot the back of the head
                plot(hn_x,hn_y,'Color',[.2 .2 .2],'LineWidth',1.2);
                plot(-hn_x,hn_y,'Color',[.2 .2 .2],'LineWidth',1.2);

                axis equal;

                % plot the channel labels
                if isfield(cfg,'plotchannellab')
                    switch cfg.plotchannellab
                        case 'yes'
                            for i = 1 : length(chanlabelb)
                                text(-xyz_b(i,2)+.01,xyz_b(i,3)+.075,chanlabelb{i},'FontSize',8);
                            end
                        otherwise
                    end
                end

                % little hack to remove some of the white space
                hold off;

                axis square off;
                set(gca,'XLim',[-1.2 1.2],'YLim',[-1.2 1.2]);

                sub_pos = get(gca,'position');
                set(gca,'Position',[sub_pos(1)-0.1 sub_pos(2)-0 sub_pos(3)*1.4 sub_pos(4)*1.4])
                text(0,-1.4,sprintf('%.0f to %.0f',cfg.latency(itopo,:)),'HorizontalAlignment','center','FontSize', 15);

            case 'top'
                % top view is where z > 0
                chanidx = find(xyz(:,3)>=-0.01);
                xyz_b = xyz(chanidx,:);
                x = -xyz_b(:,2);
                y = xyz_b(:,1);
                [xq,yq] = meshgrid(-1:meshgrn:1, -1:meshgrn:1);
                vq = griddata(x,y,v(chanidx),xq,yq,'v4');
                chanlabelb = chanlabel(chanidx);

                % trim outside the unit circle and the minimum y val
                vq((xq.^2 + yq.^2) > 1) = NaN;

                % trim outside the unit circle and the minimum y val
                % vq(yq < min(y)) = NaN;
                vq((xq.^2 + yq.^2) > 1) = NaN;

                contourf(xq,yq,vq,cfg.ncontours,'fill','on','LineStyle','none');
                if isempty(cfg.clim)
                    maxval = ceil(max(abs(v)) * 10) / 10;
                    caxis([-maxval maxval]);
                else
                    caxis(cfg.clim(:)); % add support for different clim per bin?
                end
                %colorbar;
                hold on;
                % scatter(-xyz_b(:,2),xyz_b(:,3),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
                scatter(-xyz_b(:,2),xyz_b(:,1),cfg.elecsize,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');

                if isfield(cfg,'plotchannellab')
                    switch cfg.plotchannellab
                        case 'yes'

                            for i = 1 : length(chanlabelb)
                                text(-xyz_b(i,2)+.01,xyz_b(i,1)+.075,chanlabelb{i},'FontSize',8);
                            end

                        otherwise
                    end
                end
                % head
                hn_x = -1:.01:1;
                hn_y = sqrt(1-hn_x.^2);

                % nose
                ns_x1 = -.14:grn:0;
                ns_x2 = 0:grn:.14;
                ns_y1 = 1.5*ns_x1 + 1.2;
                ns_y2 = -1.5*ns_x2 + 1.2;

                plot(hn_x, hn_y,'Color',[.2 .2 .2],'LineWidth',1.2);
                plot(hn_x,-hn_y,'Color',[.2 .2 .2],'LineWidth',1.2);
                plot(ns_x1, ns_y1,'Color',[.2 .2 .2],'LineWidth',1.2);
                plot(ns_x2, ns_y2,'Color',[.2 .2 .2],'LineWidth',1.2);
                axis equal;

                hold off;
                set(gca,'XLim',[-1.2 1.2],'YLim',[-1.2 1.2]);
                axis square off;
                sub_pos = get(gca,'position');

                % this helps to utilize more of the  space for each plot
                set(gca,'Position',[sub_pos(1)-0.1 sub_pos(2)-0 sub_pos(3)*1.4 sub_pos(4)*1.4])
                text(0,-1.2,sprintf('%.0f to %.0f',cfg.latency(itopo,:)),'HorizontalAlignment','center','FontSize', 15);


            case 'left'
                chanidx = find(xyz(:,2)>=-0.01);
                xyz_b = xyz(chanidx,:);
                chanlabelb = chanlabel(chanidx);



                x = -xyz_b(:,1);
                y = xyz_b(:,3);
                [xq,yq] = meshgrid(-1:meshgrn:1, -1:meshgrn:1);

                vq = griddata(x,y,v(chanidx),xq,yq,'v4');

                % trim outside the unit circle and the minimum y val
                % vq(yq < min(y)) = NaN;
                if isfield(cfg,'cutoff')
                    if isnumeric(cfg.cutoff)
                        cutoff = cfg.cutoff;
                    else
                        cutoff = ismember(chanlabel(chanidx),cfg.cutoff);
                        cutoff = y(cutoff);
                    end
                    vq(yq < cutoff) = NaN;
                else
                    vq((.15 .* xq.^2 + yq + .6) < 0) = NaN;
                end
                vq((xq.^2 + yq.^2) > 1) = NaN;
                vq((.15 .* xq.^2 + yq + .6) < 0) = NaN;
                vq((.5 .* (xq + 1).^2 + 0.1 + yq) < 0) = NaN;
                contourf(xq,yq,vq,cfg.ncontours,'fill','on','LineStyle','none');
                hold on;

                if isempty(cfg.clim)
                    maxval = ceil(max(abs(v)) * 10) / 10;
                    caxis([-maxval maxval]);
                else
                    caxis(cfg.clim(:)); % add support for different clim per bin?
                end
                %colorbar;
                hold on;
                % scatter(-xyz_b(:,2),xyz_b(:,3),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
                scatter(-xyz_b(:,1),xyz_b(:,3),cfg.elecsize,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');

                if isfield(cfg,'plotchannellab')
                    switch cfg.plotchannellab
                        case 'yes'

                            for i = 1 : length(chanlabelb)
                                text(-xyz_b(i,1)+.01,xyz_b(i,3)+.075,chanlabelb{i},'FontSize',8);
                            end

                        otherwise
                    end
                end

                % head & neck (hn_)
                hn_k = 0.45;     % head-neck joint
                hn_b = 1.15;     % mid-point of neck
                hn_x1 = -1:grn:hn_k;
                hn_x2 = hn_k:grn:1.2;
                hn_y1 = sqrt(1-hn_x1.^2);
                hn_a = (-hn_k / sqrt(1-hn_k.^2)) / (2*(hn_k-hn_b));     % curve of parabola
                hn_c = sqrt(1-hn_k.^2) - hn_a .* (hn_k - hn_b) .^ 2;
                hn_y2 = hn_a * (hn_x2 - hn_b) .^ 2 + hn_c;

                plot([hn_y1, hn_y2],-[hn_x1, hn_x2],'Color',[.2 .2 .2],'LineWidth',1.2);
                plot(-hn_y1(1:ceil(length(hn_x1)*.85)),-hn_x1(1:ceil(length(hn_x1)*.85)),'Color',[.2 .2 .2],'LineWidth',1.2);

                ns_y1 = -hn_x1(ceil(length(hn_x1)*.85)):-grn:-0.7;
                ns_x1 = linspace(-hn_y1(ceil(length(hn_x1)*.85)),-1.15,length(ns_y1));
                plot(ns_x1,ns_y1,'Color',[.2 .2 .2],'LineWidth',1.2);

                ns_x2 = ns_x1(end):grn:-0.9;
                ns_y2 = ns_y1(end) * ones(length(ns_x2),1);
                plot(ns_x2,ns_y2,'Color',[.2 .2 .2],'LineWidth',1.2);

                mth_x1 = -0.75:-grn:ns_x2(end);
                k = 0.5 / (mth_x1(end)-mth_x1(1)).^2;
                mth_y1 = k * (mth_x1 - mth_x1(1)).^2 - 1.2;
                plot(mth_x1,mth_y1,'Color',[.2 .2 .2],'LineWidth',1.2);

                axis equal;

                hold off;
                set(gca,'XLim',[-1.2 1.2],'YLim',[-1.2 1.2]);
                axis square off;
                sub_pos = get(gca,'position');

                % this helps to utilize more of the  space for each plot
                set(gca,'Position',[sub_pos(1)-0.1 sub_pos(2)-0 sub_pos(3)*1.4 sub_pos(4)*1.4])
                text(0,-1.2,sprintf('%.2f to %.2f',cfg.latency(itopo,:)),'HorizontalAlignment','center','FontSize', 15);

            case 'right'
                chanidx = find(xyz(:,2)<=0.01);
                xyz_b = xyz(chanidx,:);
                chanlabelb = chanlabel(chanidx);


                x = xyz_b(:,1);
                y = xyz_b(:,3);
                [xq,yq] = meshgrid(-1:meshgrn:1, -1:meshgrn:1);

                vq = griddata(x,y,v(chanidx),xq,yq,'v4');

                % trim outside the unit circle and the minimum y val
                % vq(yq < min(y)) = NaN;
                if isfield(cfg,'cutoff')
                     if isnumeric(cfg.cutoff)
                        cutoff = cfg.cutoff;
                     else
                        cutoff = ismember(chanlabel(chanidx),cfg.cutoff);
                        cutoff = y(cutoff);
                    end
                    vq(yq < cutoff) = NaN;
                else
                    vq((.15 .* xq.^2 + yq + .6) < 0) = NaN;
                end
                vq((xq.^2 + yq.^2) > 1) = NaN;
                vq((.15 .* xq.^2 + yq + .6) < 0) = NaN;
                vq((.5 .* (xq - 1).^2 + 0.1 + yq) < 0) = NaN;

                contourf(xq,yq,vq,cfg.ncontours,'fill','on','LineStyle','none');
                hold on;

                if isempty(cfg.clim)
                    maxval = ceil(max(abs(v)) * 10) / 10;
                    caxis([-maxval maxval]);
                else
                    caxis(cfg.clim(:)); % add support for different clim per bin?
                end
                %colorbar;
                hold on;
                % scatter(-xyz_b(:,2),xyz_b(:,3),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
                scatter(xyz_b(:,1),xyz_b(:,3),cfg.elecsize,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
                if isfield(cfg,'plotchannellab')
                    switch cfg.plotchannellab
                        case 'yes'

                            for i = 1 : length(chanlabelb)
                                text(xyz_b(i,1)+.01,xyz_b(i,3)+.075,chanlabelb{i},'FontSize',8);
                            end

                        otherwise
                    end
                end
                % head & neck (hn_)
                hn_k = 0.45;     % head-neck joint
                hn_b = 1.15;     % mid-point of neck
                hn_x1 = -1:grn:hn_k;
                hn_x2 = hn_k:grn:1.2;
                hn_y1 = sqrt(1-hn_x1.^2);
                hn_a = (-hn_k / sqrt(1-hn_k.^2)) / (2*(hn_k-hn_b));     % curve of parabola
                hn_c = sqrt(1-hn_k.^2) - hn_a .* (hn_k - hn_b) .^ 2;
                hn_y2 = hn_a * (hn_x2 - hn_b) .^ 2 + hn_c;

                plot(hn_y1(1:ceil(length(hn_x1)*.85)),-hn_x1(1:ceil(length(hn_x1)*.85)),'Color',[.2 .2 .2],'LineWidth',1.2);
                plot(-[hn_y1, hn_y2],-[hn_x1, hn_x2],'Color',[.2 .2 .2],'LineWidth',1.2);

                ns_y1 = -hn_x1(ceil(length(hn_x1)*.85)):-grn:-0.7;
                ns_x1 = linspace(hn_y1(ceil(length(hn_x1)*.85)),1.15,length(ns_y1));
                plot(ns_x1,ns_y1,'Color',[.2 .2 .2],'LineWidth',1.2);

                ns_x2 = .9:grn:ns_x1(end);
                ns_y2 = ns_y1(end) * ones(length(ns_x2),1);
                plot(ns_x2,ns_y2,'Color',[.2 .2 .2],'LineWidth',1.2);

                mth_x1 = 0.75:grn:ns_x2(1);
                k = 0.5 / (mth_x1(end)-mth_x1(1)).^2;
                mth_y1 = k * (mth_x1 - mth_x1(1)).^2 - 1.2;
                plot(mth_x1,mth_y1,'Color',[.2 .2 .2],'LineWidth',1.2);

                axis equal;

                hold off;
                set(gca,'XLim',[-1.2 1.2],'YLim',[-1.2 1.2]);
                axis square off;
                sub_pos = get(gca,'position');

                % this helps to utilize more of the  space for each plot
                set(gca,'Position',[sub_pos(1)-0.1 sub_pos(2)-0 sub_pos(3)*1.4 sub_pos(4)*1.4])
                text(0,-1.2,sprintf('%.0f to %.0f',cfg.latency(itopo,:)),'HorizontalAlignment','center','FontSize', 15);
        end
    end
end
