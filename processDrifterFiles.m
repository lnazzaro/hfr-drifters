function Data = processDrifterFiles(inputFile,varargin)

% input:
% inputFile - netcdf file with raw lon, lat, time
%
% output:
% Data - data structure with processed time, lon, lat, u, v, settings used
%       in processing
%
% varargin options:
% timeInterval - minutes (default: 60 minutes)
% qcRemove - array with flags to be excluded (0:no_qc_performed 1:good_data
%   2:probably_good_data 3:bad_data_that_are_potentially_correctable 
%   4:bad_data 5:value_changed 6:not_used 7:not_used 8:interpolated_value 
%   9:missing_value); empty brackets ([]) to include all values (default:
%   [3,4])
% ctdIncluded - logical indicating whether temp, cond, sal data is part of
%   dataset (default: false)
% outputDir - directory holding empty nc file to write processed data to;
%   must already have an empty file with the same file name as inputFile
%   with -processed appended to end of file name (not written to netcdf if
%   left empty)
% maxSpeed - maximum speed (cm/s) considered valid; anything exceeding this
%   in original drifter data will be excluded from processing (default: 300
%   cm/s)
% minSpeed - minimum speed (cm/s) considered valid; anything below this
%   in original drifter data will be excluded from processing (default: eps)
% maxGap - maximum time gap (hours) between points (default: 12 hours)

app = mfilename;

[~,fName,~]=fileparts(inputFile);

error_cutoff_max=300;
error_cutoff_min=eps;
outputDir='';
toNc=false;
timeInterval=60;
qcRemove=[3,4];
ctdAvail=false;
maxGap=12;

for x = 1:2:length(varargin)
    name = varargin{x};
    value = varargin{x+1};
    
    switch lower(name)
        case 'maxspeed'
            if ~isnumeric(value)&numel(value)~=1
                fprintf(2,...
                    '%s: Value for option %s must be a single numeric.\n',...
                    app,...
                    name);
                return;
            end
            error_cutoff_max = value;
        case 'minspeed'
            if ~isnumeric(value)&numel(value)~=1
                fprintf(2,...
                    '%s: Value for option %s must be a single numeric.\n',...
                    app,...
                    name);
                return;
            end
            error_cutoff_min = value;
        case 'outputdir'
            if ~isempty(value)&(~ischar(value)|~exist(value,'dir'))
                warning(...
                    '%s: Value for option %s must be an existing directory. Continuing without netcdf output.\n',...
                    app,...
                    name);
                value='';
            end
            if ~isempty(value)&~exist(fullfile(value, [fName '-processed.nc']),'file')
                warning(...
                    '%s: Empty output file %s does not exist in directory %s. Continuing without netcdf output.\n',...
                    app,...
                    [fName '-processed.nc'],...
                    value);
                value='';
            end
            if ~isempty(value)
                outputDir = value;
                toNc=true;
            end
        case 'timeinterval'
            if ~isnumeric(value)&numel(value)~=1
                fprintf(2,...
                    '%s: Value for option %s must be a single numeric.\n',...
                    app,...
                    name);
                return;
            end
            timeInterval = value;
        case 'qcremove'
            if ~isnumeric(value)
                fprintf(2,...
                    '%s: Value for option %s must be numeric.\n',...
                    app,...
                    name);
                return;
            end
            qcRemove = value;
        case 'ctdincluded'
            if ~islogical(value)|numel(value)~=1
                fprintf(2,...
                    '%s: Value for option %s must be a single logical.\n',...
                    app,...
                    name);
                return;
            end
            ctdAvail = value;
        case 'maxgap'
            if ~isnumeric(value)|numel(value)~=1
                fprintf(2,...
                    '%s: Value for option %s must be a single numeric.\n',...
                    app,...
                    name);
                return;
            end
            maxGap = value;
    end
end

flags={0,'no_qc_performed';
    1,'good_data';
    2,'probably_good_data';
    3,'bad_data_that_are_potentially_correctable';
    4,'bad_data';
    5,'value_changed';
    6,'not_used';
    7,'not_used';
    8,'interpolated_value';
    9,'missing_value'};
if(isempty(qcRemove)|isnan(qcRemove))
    removaltext='';
else
    removaltext=[', data points flagged with ' flags{qcRemove(1)+1,2}];
    for n=2:length(qcRemove)
        if(length(qcRemove)>2)
            removaltext=[removaltext ','];
        end
        if(n==length(qcRemove))
            removaltext=[removaltext ' or'];
        end
        removaltext=[removaltext ' ' flags{qcRemove(n)+1,2}];
    end
    removaltext=[removaltext ' removed before interpolation'];
end

if(toNc)
    outFile=fullfile(outputDir, [fName '-processed.nc']);
    ncwriteatt(outFile,'time','observation_type','interpolated');
    ncwriteatt(outFile,'lat','observation_type','interpolated');
    ncwriteatt(outFile,'lon','observation_type','interpolated');
    ncwriteatt(outFile,'pressure','observation_type','interpolated');
    ncwriteatt(outFile,'depth','observation_type','interpolated');
    if(ctdAvail)
        ncwriteatt(outFile,'temperature','observation_type','interpolated');
        ncwriteatt(outFile,'conductivity','observation_type','interpolated');
        ncwriteatt(outFile,'salinity','observation_type','interpolated');
        ncwriteatt(outFile,'density','observation_type','interpolated');
    end
    ncwriteatt(outFile,'u','observation_type','interpolated');
    ncwriteatt(outFile,'v','observation_type','interpolated');

    ncwriteatt(outFile,'/','date_created',...
        datestr(datetime('now')-tzoffset(datetime('now','timezone','local'))))
    atttemp=ncreadatt(inputFile,'/','processing_level');
    ncwriteatt(outFile,'/','processing_level',...
        [atttemp removaltext ', interpolated to ' int2str(timeInterval) '-minute timesteps'])
    ncwriteatt(outFile,'/','history',...
        [removaltext 'data interpolated to ' int2str(timeInterval) '-minute timesteps'])
    atttemp=ncreadatt(inputFile,'/','id');
    ncwriteatt(outFile,'/','id',atttemp);
    atttemp=ncreadatt(inputFile,'/','title');
    ncwriteatt(outFile,'/','title',[atttemp '-processed']);
    atttemp=ncread(inputFile,'trajectory');
    ncwrite(outFile,'trajectory',atttemp);
    ncwriteatt(outFile,'/','raw_data_file',[fName '.nc'])
end

% read in time, lon, lat, pressure, depth [temp, cond, sal, dens]
% read _qc from original, remove bad qc
time_orig_all=datenum(1970,1,1,0,0,double(ncread(inputFile,'time')));
time_orig=unique(time_orig_all);
time_orig=time_orig(~isnan(time_orig));
vars={'lat','lon','pressure','depth'};
if(ctdAvail)
    vars=[vars,{'temperature','conductivity','salinity','density'}];
end
for x=1:length(vars)
    v=ncread(inputFile,vars{x});
    v_qc=ncread(inputFile,[vars{x} '_qc']);
    v(ismember(v_qc,qcRemove))=nan;
    v2=nan(size(time_orig));
    for n=1:length(time_orig)
        ind=find(time_orig_all==time_orig(n));
        v2(n)=nanmedian(v(ind));
    end
    eval([vars{x} '=v2;'])
end

timeInterval=timeInterval/60/24;
time_interp=round(min(time_orig)/timeInterval)*timeInterval:timeInterval:round(max(time_orig)/timeInterval)*timeInterval;

% remove LL nans
ind=find(~isnan(lon)&~isnan(lat));
dtime=time_orig(ind);
lon=lon(ind);
lat=lat(ind);

% get distance between measurements
DI = distance(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end));
DI = deg2km(DI)*100000;
% get time difference between measurements
T = diff(dtime)*24*60*60;
% remove any over max threshold
ind=find(DI./T>error_cutoff_max|DI./T<error_cutoff_min|T>maxGap*60*60);
ind = union(ind, ind+1);  %Extraneous point crude bugfix 5/25/2007
lon(ind)=NaN;%[];
lat(ind)=NaN;%[];
dtime(ind)=NaN;%[];

% get distance on cleaned dataset
DI = distance(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end));
DI = deg2km(DI)*1000;

% get bearing (radians)
AZ =  azimuth(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end));
AZ = (90-AZ)*pi/180;

% get original time interval (s)
T = diff(dtime)*24*60*60;

% get midpoints at original times
Data.lon = (lon(1:end-1) + lon(2:end))/2;
Data.lat = (lat(1:end-1) + lat(2:end))/2;
Data.u = cos(AZ).*DI./T;
Data.v = sin(AZ).*DI./T;
Data.time_orig = (dtime(1:end-1) + dtime(2:end))/2;
Data.interval = (dtime(2:end) - dtime(1:end-1));

% interpolate
ind=find(~isnan(Data.u));
Data.time=time_interp;
Data.lon=interp1(Data.time_orig(ind),Data.lon(ind),Data.time);
Data.lat=interp1(Data.time_orig(ind),Data.lat(ind),Data.time);
Data.u=interp1(Data.time_orig(ind),Data.u(ind),Data.time);
Data.v=interp1(Data.time_orig(ind),Data.v(ind),Data.time);
diffT=diff(Data.time_orig(ind));
ind2=find(diffT>timeInterval*2);
for i=1:length(ind2)
    x=Data.time>Data.time_orig(ind(ind2(i)))+timeInterval/2&...
        Data.time<Data.time_orig(ind(ind2(i)+1))-timeInterval/2;
    Data.lon(x)=nan;
    Data.lat(x)=nan;
    Data.u(x)=nan;
    Data.v(x)=nan;
end

% interpolate the rest
for n=3:length(vars)
    dtime=time_orig;
    eval(['v=' vars{n} ';'])
    if(~all(isnan(v)))
        ind=find(~isnan(v));
        dtime=dtime(ind);
        v=v(ind);
        dtime=(dtime(1:end-1) + dtime(2:end))/2;
        v=(v(1:end-1)+v(2:end))/2;
        Data.(vars{n})=interp1(dtime,v,Data.time);
        diffT=diff(dtime);
        ind2=find(diffT>timeInterval*3/60/24);
        for i=1:length(ind2)
            x=Data.time>Data.time_orig(ind(ind2(i)))+timeInterval/2&...
                Data.time<Data.time_orig(ind(ind2(i)+1))-timeInterval/2;
            Data.(vars{n})(x)=nan;
        end
    else
        Data.(vars{n})=nan(size(Data.time));
    end
end
Data.vars=vars;

Data.attributes.rawDataFile=inputFile;
Data.attributes.max_speed=error_cutoff_max;
Data.attributes.min_speed=error_cutoff_min;
Data.attributes.ncFile=nan;
Data.attributes.timeIntervalMinutes=timeInterval*24*60;
Data.attributes.qcFlagsRemoved=qcRemove;

% write to nc
if(toNc)
    Data.attributes.ncFile=outFile;
    ncwrite(outFile,'time',(Data.time-datenum(1970,1,1))*24*60*60);
    ncwrite(outFile,'time_qc',8*ones(size(Data.time)));

    ncwrite(outFile,'u',Data.u);
    ncwrite(outFile,'u_qc',8*ones(size(Data.u)));
    ncwrite(outFile,'v',Data.v);
    ncwrite(outFile,'v_qc',8*ones(size(Data.v)));

    for n=1:length(vars)
        v=vars{n};
        ncwrite(outFile,vars{n},Data.(vars{n}));
        if(~all(isnan(Data.(vars{n}))))
            ncwrite(outFile,[vars{n} '_qc'],8*ones(size(Data.time)));
        end
    end
end

Data.u=Data.u*100;
Data.v=Data.v*100;
Data.units='cm/s';