function stats=drifterVradials(radialStructs,varargin)

% input:
% radialStructs - single structured array or cell array containing multiple
%       structured arrays of HFR radial velocities and matching drifter
%       velocities (rotated), as output from one or multiple iterations of 
%       drifter2hfr.m; each should be only one dataset 
%       (ie radialSpeeds.LOVE.measured, etc)
%
% output:
% stats - structured array with statistics for each dataset, including
%       number of comparison points (N), RMSE, and correlation
%       coefficient (r)
% images include map of drifter track (if drifter variable provided as
%       input) and line plot with all radial velocities and matching
%       drifter velocities in lighter shades
%
% varargin options:
% plot - logical indicating whether or not to generate plots (default:
%       true)
% names - cell array of names for each radial dataset input (defaults to
%       #1, #2, etc)
% colors - cell array of colors to use in line plots for each radial 
%       dataset (defaults to blue, red, black, green), in order, can be
%       matlab color strings or rgb values
% t0 - start time
% t1 - end time
% drifter - structured array with drifter data; only needed if you want to
%       plot the drifter track
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

app = mfilename;

if ~iscell(radialStructs)
    radialStructs={radialStructs};
end

stats(length(radialStructs))=struct('name',[],'N',[],'RMSE',[],'r',[]);

plotFigs=true;
colors={[0 0 1],[1 0 0],[0 0 0],[0 1 0]}; % blue, red, black, green
colors=colors(1:length(radialStructs));
names=cell(length(radialStructs),1);
for n=1:length(names)
    names{n}=num2str(n);
end
ind=~isnan(radialStructs{1}.rotated_drifter_velocity) & ...
    ~isnan(radialStructs{1}.HFR_radial_velocity);
t0=min(radialStructs{1}.time(ind));
t1=max(radialStructs{1}.time(ind));
for n=2:length(radialStructs)
    ind=~isnan(radialStructs{n}.rotated_drifter_velocity) & ...
        ~isnan(radialStructs{n}.HFR_radial_velocity);
    if(sum(ind)>0)
        t0=nanmin(t0,min(radialStructs{n}.time(ind)));
        t1=nanmax(t1,max(radialStructs{n}.time(ind)));
    end
end
drifter=struct([]);
bathyDir=[];
bathyFile=[];
isobaths=[-20,-50,-100];

if(isempty(t0)|isempty(t1))
    fprintf(2,...
        '%s: No non-NaN data in file.\n',...
        app);
    return;
end
    

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
            if ~iscell(value) | length(value)~=length(radialStructs) | ~all(cellfun(@ischar,value))
                fprintf(2,...
                    '%s: Value for option %s must be a cell array the same length as radialStructs containing only strings.\n',...
                    app,...
                    name);
                return;
            end
            names=value;
        case 'colors'
            if ~iscell(value) | length(value)~=length(radialStructs)
                fprintf(2,...
                    '%s: Value for option %s must be a cell array the same length as radialStructs.\n',...
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
        case 'drifter'
            if ~isstruct(value)
                fprintf(2,...
                    '%s: Value for option %s must be a struct.\n',...
                    app,...
                    name);
                return;
            end
            drifter=value;
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


for n=1:length(radialStructs)
    stats(n).name=names{n};
    ind=~isnan(radialStructs{n}.rotated_drifter_velocity) & ...
        ~isnan(radialStructs{n}.HFR_radial_velocity) & ...
        radialStructs{n}.time>=t0 & radialStructs{n}.time<=t1;
    stats(n).N=sum(ind);
    stats(n).RMSE=sqrt(sum((radialStructs{n}.HFR_radial_velocity(ind)-...
        radialStructs{n}.rotated_drifter_velocity(ind)).^2)/sum(ind));
    c=corrcoef(radialStructs{n}.rotated_drifter_velocity(ind),...
        radialStructs{n}.HFR_radial_velocity(ind));
    stats(n).r=c(2);
end


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

if ~isempty(drifter)
    figure
    ind=find(drifter.time>=t0&drifter.time<=t1&...
        ~isnan(drifter.lon)&~isnan(drifter.lat));
    xl=[min(drifter.lon(ind)) max(drifter.lon(ind))]+[-2 2];
    yl=[min(drifter.lat(ind)) max(drifter.lat(ind))]+[-2 2];
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
    plot(drifter.lon(ind),drifter.lat(ind),'k');
    scatter(drifter.lon(ind),drifter.lat(ind),5,drifter.time(ind),'filled');
    grid on
end

figure
hold on

legnames={};

for n=1:length(radialStructs)
    c_hfr=colors{n};
    cd=1-c_hfr;
    c_drifter=c_hfr+.75*cd;
    plot(radialStructs{n}.time,radialStructs{n}.rotated_drifter_velocity,'col',c_drifter,'marker','.','linewidth',2)
    plot(radialStructs{n}.time,radialStructs{n}.HFR_radial_velocity,'col',c_hfr,'marker','.')
    legnames=[legnames,{['Drifter ' names{n}]},{['HFR ' names{n}]}];
end

set(gca,'xtick',t0tick:tint:t1tick,'xticklabel',datestr(t0tick:tint:t1tick,tformat))
grid on
xlim([t0 t1])
legend(legnames)

