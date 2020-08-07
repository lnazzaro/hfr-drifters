function stats=drifterVtotals(drifterStruct,totalsStructs,varargin)

app = mfilename;

if ~iscell(totalsStructs)
    totalsStructs={totalsStructs};
end

stats(length(totalsStructs))=struct('name',[],'N',[],'RMSEu',[],'RMSEv',[],'corr_mag',[],'corr_dir',[]);

plotFigs=true;
% colors={[0 0 1],[1 0 0],[0 0 0],[0 1 0]}; % blue, red, black, green
% colors=colors(1:length(totalsStructs));
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
%         case 'colors'
%             if ~iscell(value) | length(value)~=length(totalsStructs)
%                 fprintf(2,...
%                     '%s: Value for option %s must be a cell array the same length as totalsStructs.\n',...
%                     app,...
%                     name);
%                 return;
%             end
%             for n=1:length(value)
%                 if ismember(value{n},{'y','yellow'})
%                     value{n}=[1 1 0];
%                 elseif ismember(value{n},{'m','magenta'})
%                     value{n}=[1 0 1];
%                 elseif ismember(value{n},{'c','cyan'})
%                     value{n}=[0 1 1];
%                 elseif ismember(value{n},{'r','red'})
%                     value{n}=[1 0 0];
%                 elseif ismember(value{n},{'g','green'})
%                     value{n}=[0 1 0];
%                 elseif ismember(value{n},{'b','blue'})
%                     value{n}=[0 0 1];
%                 elseif ismember(value{n},{'w','white'})
%                     value{n}=[1 1 1];
%                 elseif ismember(value{n},{'k','black'})
%                     value{n}=[0 0 0];
%                 elseif ~isnumeric(value{n})|length(value{n})~=3|min(value{n})<0|max(value{n})>1
%                     fprintf(2,...
%                         '%s: Value for option %s is not a recognized color or rgb value (must range 0-1).\n',...
%                         app,...
%                         name);
%                     return;
%                 end
%             end
%             colors=value;
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


if ~plotFigs
    return
end


tformat='mm/dd';
if t1-t0<3
    tint=.5;
    tformat='mm/dd HH';
elseif t1-t0<8
    tint=1
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
caxis([t0 t1])
set(cb,'ytick',t0tick:tint:t1tick,'yticklabel',datestr(t0tick:tint:t1tick,tformat))
plot(drifterStruct.lon(ind),drifterStruct.lat(ind),'k');
scatter(drifterStruct.lon(ind),drifterStruct.lat(ind),5,drifterStruct.time(ind),'filled');
grid on


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

end


% feather plots

nFigs=ceil(t1-t0)/7;
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
    
    subplot(length(totalsStructs)+1,1,1)
    ind=find(drifterStruct.time>=t0f&drifterStruct.time<=t1f);
    feather(drifterStruct.u(ind),drifterStruct.v(ind),'k')
    hold on
    plot([1 length(ind)],[0 0],'color',[.5 .5 .5])
    iticks=ismember(drifterStruct.time(ind),tticks);
    iticklabels=datestr(drifterStruct.time(ind(iticks)),tformat);
    set(gca,'xtick',find(iticks),'xticklabel',iticklabels,'xminortick','on')
    grid(gca,'minor')
    grid on
    ylabel('Velocity (cm/s)')
    xlabel(['Time (' tformat ')'])
    xlim([-2 length(ind)+3])
    title(['Drifter ' datestr(t0f) '-' datestr(t1f)])
    
    for n=1:length(totalsStructs)
        subplot(length(totalsStructs)+1,1,n+1)
        ind=find(totalsStructs{n}.time>=t0f&totalsStructs{n}.time<=t1f);
        feather(totalsStructs{n}.HFR_totals_u(ind),totalsStructs{n}.HFR_totals_v(ind),'k')
        hold on
        plot([1 length(ind)],[0 0],'color',[.5 .5 .5])
        iticks=ismember(totalsStructs{n}.time(ind),tticks);
        iticklabels=datestr(totalsStructs{n}.time(ind(iticks)),tformat);
        set(gca,'xtick',find(iticks),'xticklabel',iticklabels,'xminortick','on')
        grid(gca,'minor')
        grid on
        ylabel('Velocity (cm/s)')
        xlabel(['Time (' tformat ')'])
        xlim([-2 length(ind)+3])
        title(names{n})
    end
    
end
    