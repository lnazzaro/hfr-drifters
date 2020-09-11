function [stats,figureinfo]=drifterVtotals(drifterStruct,totalsStructs,varargin)

% input:
% totalsStructs - single structured array or cell array containing multiple
%       structured arrays of HFR totals velocities, as output from one or 
%       multiple iterations of drifter2hfr.m; each should be only one 
%       dataset (ie totalsVelocities.totals, etc)
%
% output:
% stats - structured array with statistics for each dataset, including
%       number of comparison points (N), RMSE of u velocities (RMSEu), 
%       RMSE of v velocities (RMSEv), complex correlation strength
%       (corr_mag), and complex correlation directional offset (corr_dir)
% images include map of drifter track, "current roses" showing direction
%       TOWARDS over given time period for the drifter and the matching HFR
%       and directional offset between drifter and HFR, and subsets of the
%       drifter track with drifter and HFR velocities at each timestamp and
%       location, for each totals dataset
%
% varargin options:
% plot - logical indicating whether or not to generate plots (default:
%       true)
% names - cell array of names for each totals dataset input (defaults to
%       #1, #2, etc)
% colors - cell array of colors to use in quiver plots for each totals 
%       dataset (defaults to blue, red, green), in order, can be
%       matlab color strings or rgb values; drifter velocity will be black
% t0 - start time
% t1 - end time
% bathymetryFile - netcdf file with lon, lat, altitude (only used for
%       drifter track plot)
% bathymetryDir - directory containing bathymetry netcdf files holding lon,
%       lat, altitude; file names should contain lower left and upper right
%       domain bounds and resolution, as in 
%       GMRTv3_6_20190703_LL_76W_38N_UR_72W_42N_res_244m.grd
%       (ignored if bathymetryFile is provided, only used for drifter track
%       plot)
% isobaths - isobaths to plot with drifter track (only used for
%       drifter track plot)
% quiverScaleFactor - scale factor for plotting drifter velocities on map
%       (default 0.01)
% plotSpeedDifference - plot all velocities at all locations with speed
%       difference between drifter and HFR exceeding X cm/s (default 10),
%       regardless of subsetting
% plotDirectionDifference - plot all velocities at all locations with
%       direction difference between drifter and HFR exceeding X degrees 
%       (default 45), regardless of subsetting
% quiverInterval - subsetting interval for plotting velocities on map
%       (default every point)
% maxDaysInPlot - longest timespan (days) to plot in each separate drifter
%       track w/ current vectors figure


app = mfilename;

if ~iscell(totalsStructs)
    totalsStructs={totalsStructs};
end

stats(length(totalsStructs))=struct('name',[],'N',[],'RMSEu',[],'RMSEv',[],'corr_mag',[],'corr_dir',[]);

plotFigs=true;
colors={[0 0 1],[1 0 0],[0 1 0]}; % blue, red, green
colors=colors(1:length(totalsStructs));
names=cell(length(totalsStructs),1);
for n=1:length(names)
    names{n}=['HFR Totals #' num2str(n)];
end
ind=~isnan(drifterStruct.u) & ...
    ~isnan(drifterStruct.v);
t0=min(drifterStruct.time(ind));
t1=max(drifterStruct.time(ind));
for n=1:length(totalsStructs)
    ind=~isnan(totalsStructs{n}.HFR_totals_u) & ...
        ~isnan(totalsStructs{n}.HFR_totals_v);
    t0=min(t0,min(totalsStructs{n}.time(ind)));
    t1=max(t1,max(totalsStructs{n}.time(ind)));
end
bathyDir=[];
bathyFile=[];
isobaths=[-20,-50,-100];
sf=.001;
minSpdDiff=10;
minDirDiff=45;
qint=1;
maxDays=2;

for x = 1:2:length(varargin)
    name = varargin{x};
    value = varargin{x+1};
    
    switch lower(name)
        case 'plot'
            if ~islogical(value) | numel(value)~=1
                fprintf(2,...
                    '%s: Value for option %s must be a logical.\n',...
                    app,...
                    name);
                return;
            end
            plotFigs=value;
        case 'names'
            if ~iscell(value) | length(value)~=length(totalsStructs) | ~all(cellfun(@ischar,value))
                fprintf(2,...
                    '%s: Value for option %s must be a cell array the same length as totalsStructs containing only strings.\n',...
                    app,...
                    name);
                return;
            end
            names=value;
        case 'colors'
            if ~iscell(value) | length(value)~=length(totalsStructs)
                fprintf(2,...
                    '%s: Value for option %s must be a cell array the same length as totalsStructs.\n',...
                    app,...
                    name);
                return;
            end
            for n=1:length(value)
                if ismember(value{n},{'y','yellow'})
                    value{n}=[1 1 0];
                elseif ismember(value{n},{'m','magenta'})
                    value{n}=[1 0 1];
                elseif ismember(value{n},{'c','cyan'})
                    value{n}=[0 1 1];
                elseif ismember(value{n},{'r','red'})
                    value{n}=[1 0 0];
                elseif ismember(value{n},{'g','green'})
                    value{n}=[0 1 0];
                elseif ismember(value{n},{'b','blue'})
                    value{n}=[0 0 1];
                elseif ismember(value{n},{'w','white'})
                    value{n}=[1 1 1];
                elseif ismember(value{n},{'k','black'})
                    value{n}=[0 0 0];
                elseif ~isnumeric(value{n})|length(value{n})~=3|min(value{n})<0|max(value{n})>1
                    fprintf(2,...
                        '%s: Value for option %s is not a recognized color or rgb value (must range 0-1).\n',...
                        app,...
                        name);
                    return;
                end
            end
            colors=value;
        case 't0'
            if ~isnumeric(value)|value<datenum(1995,1,1)|value>now+10
                fprintf(2,...
                    '%s: Value for option %s must be a MATLAB datenum.\n',...
                    app,...
                    name);
                return;
            end
            t0=value;
        case 't1'
            if ~isnumeric(value)|value<datenum(1995,1,1)|value>now+10
                fprintf(2,...
                    '%s: Value for option %s must be a MATLAB datenum.\n',...
                    app,...
                    name);
                return;
            end
            t1=value;
        case 'bathymetrydir'
            if ~isdir(value)
                fprintf(2,...
                    '%s: Value for option %s must be a directory.\n',...
                    app,...
                    name);
                return;
            end
            bathyDir=value;
        case 'bathymetryfile'
            if ~ischar(value)
                fprintf(2,...
                    '%s: Value for option %s must be a file.\n',...
                    app,...
                    name);
                return;
            end
            bathyFile=value;
        case 'isobaths'
            if ~isnumeric(value)
                fprintf(2,...
                    '%s: Value for option %s must be numeric.\n',...
                    app,...
                    name);
                return;
            end
            isobaths=-abs(value);
        case 'quiverscalefactor'
            if ~isnumeric(value)|numel(value)>1
                fprintf(2,...
                    '%s: Value for option %s must be single-element numeric.\n',...
                    app,...
                    name);
                return;
            end
            sf=value;
        case 'plotspeeddifference'
            if ~isnumeric(value)|numel(value)>1
                fprintf(2,...
                    '%s: Value for option %s must be single-element numeric.\n',...
                    app,...
                    name);
                return;
            end
            minSpdDiff=value;
        case 'plotdirectiondifference'
            if ~isnumeric(value)|numel(value)>1
                fprintf(2,...
                    '%s: Value for option %s must be single-element numeric.\n',...
                    app,...
                    name);
                return;
            end
            minDirDiff=value;
        case 'quiverinterval'
            if ~isnumeric(value)|numel(value)>1
                fprintf(2,...
                    '%s: Value for option %s must be single-element numeric.\n',...
                    app,...
                    name);
                return;
            end
            qint=value;
        case 'maxdaysinplot'
            if ~isnumeric(value)|numel(value)>1
                fprintf(2,...
                    '%s: Value for option %s must be single-element numeric.\n',...
                    app,...
                    name);
                return;
            end
            maxDays=value;
    end
end


for n=1:length(totalsStructs)
    stats(n).name=names{n};
    indTotals=find(~isnan(totalsStructs{n}.HFR_totals_u) & ...
        ~isnan(totalsStructs{n}.HFR_totals_u) & ...
        totalsStructs{n}.time>=t0 & totalsStructs{n}.time<=t1);
    indDrifter=find(~isnan(drifterStruct.u) & ...
        ~isnan(drifterStruct.v) & ...
        drifterStruct.time>=t0 & drifterStruct.time<=t1);
    [~,indD,indHFR]=intersect(drifterStruct.time(indDrifter),totalsStructs{n}.time(indTotals));
    indD=indDrifter(indD);
    indHFR=indTotals(indHFR);
    stats(n).N=length(indD);
    stats(n).RMSEu=sqrt(sum((totalsStructs{n}.HFR_totals_u(indHFR)-...
        drifterStruct.u(indD)).^2)/length(indD));
    stats(n).RMSEv=sqrt(sum((totalsStructs{n}.HFR_totals_v(indHFR)-...
        drifterStruct.v(indD)).^2)/length(indD));
    [c, d]=complexCorr(length(indD),...
        drifterStruct.u(indD),drifterStruct.v(indD),...
        totalsStructs{n}.HFR_totals_u(indHFR),totalsStructs{n}.HFR_totals_v(indHFR));
    stats(n).corr_mag=c;
    stats(n).corr_dir=d;
end

figureinfo=[];

if ~plotFigs
    return
end


tformat='mm/dd';
if t1-t0<3
    tint=.5;
    tformat='mm/dd HH';
elseif t1-t0<8
    tint=1;
else
    tint=ceil((t1-t0)/8);
end
t0tick=ceil(t0/tint)*tint;
t1tick=floor(t1/tint)*tint;

figure
ind=find(drifterStruct.time>=t0&drifterStruct.time<=t1&...
    ~isnan(drifterStruct.lon)&~isnan(drifterStruct.lat));
xl=[min(drifterStruct.lon(ind)) max(drifterStruct.lon(ind))]+[-2 2];
yl=[min(drifterStruct.lat(ind)) max(drifterStruct.lat(ind))]+[-2 2];
if ~isempty(bathyFile) & ~isempty(bathyDir)
    bathyFile=fullfile(bathyDir,bathyFile);
elseif ~isempty(bathyDir)
    bathyFiles=dir([bathyDir '*.grd']);
    minlon=nan(1,length(bathyFiles));
    maxlon=nan(1,length(bathyFiles));
    minlat=nan(1,length(bathyFiles));
    maxlat=nan(1,length(bathyFiles));
    res=nan(1,length(bathyFiles));
    for n=1:length(bathyFiles)
        pieces=split(bathyFiles(n).name,'_');
        i=find(strcmp(pieces,'LL'));
        for p=1:2
            loc=pieces{i+p};
            if loc(end)=='W'
                minlon(n)=-str2double(loc(1:end-1));
            elseif loc(end)=='E'
                minlon(n)=str2double(loc(1:end-1));
            elseif loc(end)=='N'
                minlat(n)=str2double(loc(1:end-1));
            elseif loc(end)=='S'
                minlat(n)=-str2double(loc(1:end-1));
            end
        end
        i=find(strcmp(pieces,'UR'));
        for p=1:2
            loc=pieces{i+p};
            if loc(end)=='W'
                maxlon(n)=-str2double(loc(1:end-1));
            elseif loc(end)=='E'
                maxlon(n)=str2double(loc(1:end-1));
            elseif loc(end)=='N'
                maxlat(n)=str2double(loc(1:end-1));
            elseif loc(end)=='S'
                maxlat(n)=-str2double(loc(1:end-1));
            end
        end
        res(n)=str2double(pieces{end}(1:end-5));
    end
    bathyind=find(minlon<xl(1));
    bathyind=intersect(bathyind,find(maxlon>xl(2)));
    bathyind=intersect(bathyind,find(minlat<yl(1)));
    bathyind=intersect(bathyind,find(maxlat>yl(2)));
    if isempty(bathyind)
        bathyFile=[];
    else
        [~,bestbathy]=min(res(bathyind));
        bathyFile=[bathyDir bathyFiles(bathyind(bestbathy)).name];
    end
end
try
    bathylon=ncread(bathyFile,'lon');
    bathylat=ncread(bathyFile,'lat');
    bathyd=ncread(bathyFile,'altitude');
catch
    bathyFile=[];
end


if isempty(bathyFile)
    geoshow('landareas.shp', 'FaceColor', [.95 .88 .66]);
    hold on
    xlim(xl)
    ylim(yl)
else
    contour(bathylon,bathylat,bathyd',[0 0],'k','linewidth',1.2)
    hold on
    if ~isempty(isobaths)
        contour(bathylon,bathylat,bathyd',isobaths,'color',[.5 .5 .5])
    end
    xlim(xl)
    ylim(yl)
    project_mercator
end

cb=colorbar;
colormap(jet)
caxis([t0 t1])
set(cb,'ytick',t0tick:tint:t1tick,'yticklabel',datestr(t0tick:tint:t1tick,tformat))
plot(drifterStruct.lon(ind),drifterStruct.lat(ind),'k');
scatter(drifterStruct.lon(ind),drifterStruct.lat(ind),5,drifterStruct.time(ind),'filled');
grid on

figureinfo=[figureinfo,{'full_drifter_track'}];


% "current rose" (wind rose)
% histograms/"current rose" of separation in speed and direction

for n=1:length(totalsStructs)

    figure
    set(gcf,'position',[100 250 1200 400])
    hold on
    
    indTotals=find(~isnan(totalsStructs{n}.HFR_totals_u) & ...
        ~isnan(totalsStructs{n}.HFR_totals_u) & ...
        totalsStructs{n}.time>=t0 & totalsStructs{n}.time<=t1);
    indDrifter=find(~isnan(drifterStruct.u) & ...
        ~isnan(drifterStruct.v) & ...
        drifterStruct.time>=t0 & drifterStruct.time<=t1);
    [~,indD,indHFR]=intersect(drifterStruct.time(indDrifter),totalsStructs{n}.time(indTotals));
    indD=indDrifter(indD);
    indHFR=indTotals(indHFR);
    uD=drifterStruct.u(indD);
    vD=drifterStruct.v(indD);
    uHFR=totalsStructs{n}.HFR_totals_u(indHFR);
    vHFR=totalsStructs{n}.HFR_totals_v(indHFR);
    
    speedD=sqrt(uD.^2+vD.^2);
    dirD=atan2(uD,vD)*180/pi;
    speedHFR=sqrt(uHFR.^2+vHFR.^2);
    dirHFR=atan2(uHFR,vHFR)*180/pi;
    speedInt=ceil(max([max(speedD),max(speedHFR)])/7);
    
    subplot(1,3,1)
    HFR_rose(dirD,speedD,'dtype','meteo','di',0:speedInt:speedInt*7,'n',32,'parent',gca);
    title('Drifter')
    
    subplot(1,3,2)
    HFR_rose(dirHFR,speedHFR,'dtype','meteo','di',0:speedInt:speedInt*7,'n',32,'parent',gca);
    title(names{n})
    
    subplot(1,3,3)
    HFR_rose(dirHFR-dirD,speedHFR-speedD,'dtype','meteo','nsewlabel','meteo','n',32,'parent',gca);
    title([names{n} ' minus Drifter'])
    
    figureinfo=[figureinfo,{['current_rose_drifter_vs_' names{n}]}];

end


% point current plots

nFigs=ceil((t1-t0)/maxDays);
tInt=(t1-t0)/nFigs;

for k=1:nFigs
    
    figure
    hold on
    
    t0f=t0+tInt*(k-1);
    t1f=t0+tInt*k;
    tif=t1f-t0f;
    if tif<.5
        tticks=t0f:1/12:t1f;
        tformat='HH:MM';
    elseif tif<1
        tticks=t0f:1/8:t1f;
        tformat='HH:MM';
    elseif tif<2
        tticks=t0f:1/3:t1f;
        tformat='dd HH:MM';
    elseif tif<4
        tticks=ceil(t0f*2)/2:1/2:floor(t1f*2)/2;
        tformat='dd HH:MM';
    else 
        tticks=ceil(t0f):floor(t1f);
        tformat='dd HH:MM';
    end
    
    ind=find(drifterStruct.time>=t0f&drifterStruct.time<=t1f);
    tsub=drifterStruct.time(ind);
    drifterLon=drifterStruct.lon(ind);
    drifterLat=drifterStruct.lat(ind);
    drifterU=drifterStruct.u(ind);
    drifterV=drifterStruct.v(ind);
    drifterSpd=sqrt(drifterU.^2+drifterV.^2);
    drifterDir=atan2(drifterV,drifterU)*180/pi;
    drifterDir=mod(drifterDir,360);
    
    for n=1:length(totalsStructs)
        ind=find(totalsStructs{n}.time>=t0f&totalsStructs{n}.time<=t1f);
        totalsU{n}=totalsStructs{n}.HFR_totals_u(ind);
        totalsV{n}=totalsStructs{n}.HFR_totals_v(ind);
        totalsSpd{n}=sqrt(totalsU{n}.^2+totalsV{n}.^2);
        totalsDir{n}=atan2(totalsV{n},totalsU{n})*180/pi;
        totalsDir{n}=mod(totalsDir{n},360);
        diffSpd{n}=abs(totalsSpd{n}-drifterSpd);
        diffDir{n}=abs(totalsDir{n}-drifterDir);
        diffDir{n}(diffDir{n}>180)=360-diffDir{n}(diffDir{n}>180);
    end
    
    xl=[min(drifterLon) max(drifterLon)]+[-.1 .1];
    yl=[min(drifterLat) max(drifterLat)]+[-.1 .1];
    
    if isempty(bathyFile)
        geoshow('landareas.shp', 'FaceColor', [.95 .88 .66]);
        hold on
        xlim(xl)
        ylim(yl)
    else
        contour(bathylon,bathylat,bathyd',[0 0],'k','linewidth',1.2)
        hold on
        if ~isempty(isobaths)
            contour(bathylon,bathylat,bathyd',isobaths,'color',[.5 .5 .5])
        end
        xlim(xl)
        ylim(yl)
        project_mercator
    end

    cb=colorbar;
    colormap(jet)
    caxis([t0f t1f])
    set(cb,'ytick',tticks,'yticklabel',datestr(tticks,tformat))
    scatter(drifterLon,drifterLat,50,tsub,'filled');
    grid on
    
    
    indna=find(~isnan(drifterU)&~isnan(drifterV));
    indSpd=[];
    indDir=[];
    for n=1:length(totalsStructs)
        indSpd=union(indSpd,find(diffSpd{n}>minSpdDiff));
        indDir=union(indDir,find(diffDir{n}>minDirDiff));
    end
    indPlot=union(indSpd,indDir);
    indPlot=union(indPlot,1:qint:length(drifterLon));
    q=quiver(drifterLon(union(indPlot,indna)),drifterLat(union(indPlot,indna)),...
        drifterU(union(indPlot,indna))./cosd(drifterLat(union(indPlot,indna)))*sf,...
        drifterV(union(indPlot,indna))*sf,...
        0,'k');
    for n=1:length(totalsStructs)
        indna=find(~isnan(totalsU{n})&~isnan(totalsV{n}));
        q=[q,quiver(drifterLon(union(indPlot,indna)),...
            drifterLat(union(indPlot,indna)),...
            totalsU{n}(union(indPlot,indna))./cosd(drifterLat(union(indPlot,indna)))*sf,...
            totalsV{n}(union(indPlot,indna))*sf,...
            0,'color',colors{n})];
    end
    legend(q,['Drifter',names])
    
    title([datestr(t0f) ' - ' datestr(t1f)])
    
    figureinfo=[figureinfo,{['point_currents_' datestr(t0f,'yyyymmddTHHMM') '-' datestr(t1f,'yyyymmddTHHMM')]}];
    
end


% % feather plots
% 
% nFigs=ceil(t1-t0)/7;
% tInt=(t1-t0)/nFigs;
% 
% for k=1:nFigs
%     
%     figure
%     hold on
%     
%     t0f=t0+tInt*(k-1);
%     t1f=t0+tInt*k;
%     tif=t1f-t0f;
%     if tif<.5
%         tticks=t0f:1/12:t1f;
%         tformat='HH:MM';
%     elseif tif<1
%         tticks=t0f:1/8:t1f;
%         tformat='HH:MM';
%     elseif tif<2
%         tticks=t0f:1/3:t1f;
%         tformat='dd HH:MM';
%     elseif tif<4
%         tticks=ceil(t0f*2)/2:1/2:floor(t1f*2)/2;
%         tformat='dd HH:MM';
%     else 
%         tticks=ceil(t0f):floor(t1f);
%         tformat='dd HH:MM';
%     end
%     
%     subplot(length(totalsStructs)+1,1,1)
%     ind=find(drifterStruct.time>=t0f&drifterStruct.time<=t1f);
%     feather(drifterStruct.u(ind),drifterStruct.v(ind),'k')
%     hold on
%     plot([1 length(ind)],[0 0],'color',[.5 .5 .5])
%     iticks=ismember(round(drifterStruct.time(ind)*24),round(tticks*24));
%     iticklabels=datestr(drifterStruct.time(ind(iticks)),tformat);
%     set(gca,'xtick',find(iticks),'xticklabel',iticklabels,'xminortick','on')
%     grid(gca,'minor')
%     grid on
%     ylabel('Velocity (cm/s)')
%     xlabel(['Time (' tformat ')'])
%     xlim([-2 length(ind)+3])
%     title(['Drifter ' datestr(t0f) '-' datestr(t1f)])
%     
%     for n=1:length(totalsStructs)
%         subplot(length(totalsStructs)+1,1,n+1)
%         ind=find(totalsStructs{n}.time>=t0f&totalsStructs{n}.time<=t1f);
%         feather(totalsStructs{n}.HFR_totals_u(ind),totalsStructs{n}.HFR_totals_v(ind),'k')
%         hold on
%         plot([1 length(ind)],[0 0],'color',[.5 .5 .5])
%         iticks=ismember(round(totalsStructs{n}.time(ind)*24),round(tticks*24));
%         iticklabels=datestr(totalsStructs{n}.time(ind(iticks)),tformat);
%         set(gca,'xtick',find(iticks),'xticklabel',iticklabels,'xminortick','on')
%         grid(gca,'minor')
%         grid on
%         ylabel('Velocity (cm/s)')
%         xlabel(['Time (' tformat ')'])
%         xlim([-2 length(ind)+3])
%         title(names{n})
%     end
%     
% end
%     
