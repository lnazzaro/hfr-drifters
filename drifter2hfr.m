function [drifterVelocities, radialSpeeds, totalsVelocities] = ...
    drifter2hfr(processedDrifterData,varargin)

% input:
% processedDrifterData - netcdf file with processed time, lon, lat, u, v OR 
%       structured array output from processDrifterFiles.m, at timestamps
%       matching HFR data
%
% output:
% drifterVelocities - structured array with processed drifter data (time,
%       lon, lat, u, v) and u,v units
% radialSpeeds - structured array with information on settings used and
%       nearest HFR velocity and drifter velocity rotated based on radial
%       bearing, as well as any other associated radial data for all radial
%       datasets (as radialSpeeds.SITE.type)
% totalsVelocities  -structured array with information on settings used and
%       HFR totals u and v matched to drifter velocities as in 
%       drifterVelocities, as well as any other associated data for all 
%       totals datasets (as totalsVelocities.totals)
%
% radialSites: cell array of radial site codes to compare drifter to (or
%       'all' for all within range or 'none' (default: all)
% radialType: measured, ideal, or both (default: both)
% totalsNetworks: which totals network(s) to compare to or 'all' or 'none'
%       (only applicable if using thredds; default: all)
% radialDir: directory to pull radial data from; should have subdirectories
%       for each site code (default: /home/codaradm/data/radials/)
% totalsDir: directory with totals files to compare to, or 'thredds'
% saveAllFields: logical (default true) indicating whether to save only u
%       and v, or all variables available in comparison dataset
% outFile: name of netcdf file (including directory) to save output to
%       (default: none); if same as processedDrifterFile, new fields will
%       be added to existing file - STILL UNDER DEVELOPMENT
% maxSeparation: maximum distance (km) from drifter to consider in
%       comparison (default: 10)
% totalsComparisonMethod: how to match to totals data; options include
%       nearest, mean, median, weightedgaussian, weightedexponential,
%       weightedgaussianellipse, and weightedexponentialellipse (default:
%       weightedgaussian); ellipses are based on relative u and v velocity
%       calculated from drifter data
% totalsToRemove: QC rules indicating totals vectors to remove as cell 
%       array (default: u_err>0.6 or v_err>0.6, given as
%       {'u_err>0.6','v_err>0.6'}
% radialsToRemove: QC rules indicating radial vectors to remove as cell 
%       array (default: {'VFLG>0'}); must use 4-character column names


app = mfilename;


radialSites={'all'};
totalsNetworks={'all'};
radialDir='/home/codaradm/data/radials/';
totalsDir='thredds';
radialType='both';
%totalsDir='/home/codaradm/data/totals/maracoos/oi/nc/'; % 5MHz/
saveAllFields=true;
outFile='none';
maxSeparation=10;
totalsComparisonMethod='weightedgaussian';
totalsToRemove={'u_err>0.6','v_err>0.6'};
radialsToRemove={'VFLG>0'};

totalsMethodOptions={'nearest','mean','median',...
    'weightedgaussian','weightedexponential',...
    'weightedgaussianellipse','weightedexponentialellipse'};

for x = 1:2:length(varargin)
    name = varargin{x};
    value = varargin{x+1};
    
    switch lower(name)
        case 'radialsites'
            if ~iscell(value)
                value={value};
            end
            if ~all(cellfun(@ischar,value))
                fprintf(2,...
                    '%s: Value for option %s must be a cell array containing only strings.\n',...
                    app,...
                    name);
                return;
            end
            radialSites = value;
        case 'radialtype'
            if ~ischar(value) | ~ismember(lower(value),{'ideal','measured','both','none'})
                fprintf(2,...
                    '%s: Value for option %s must be a string: ideal, measured, none, or both.\n',...
                    app,...
                    name);
                return;
            end
            radialType=lower(value);
        case 'totalsnetworks'
            if ~iscell(value)
                value={value};
            end
            if ~all(cellfun(@ischar,value))
                fprintf(2,...
                    '%s: Value for option %s must be a cell array containing only strings.\n',...
                    app,...
                    name);
                return;
            end
            totalsNetworks = value;
        case 'radialdir'
            if ~isdir(value)
                fprintf(2,...
                    '%s: Value for option %s must be a directory.\n',...
                    app,...
                    name);
                return;
            end
            if value(end)~='/'
                value=[value '/'];
            end
            radialDir = value;
        case 'totalsdir'
            if ~isdir(value) & ~strcmpi(value,'thredds')
                fprintf(2,...
                    '%s: Value for option %s must be a directory or ''thredds'' to use thredds server.\n',...
                    app,...
                    name);
                return;
            end
            if value(end)~='/' & ~strcmpi(value,'thredds')
                value=[value '/'];
            end
            totalsDir = value;
        case 'saveallfields'
            if numel(value)~=1 | ~islogical(value)
                fprintf(2,...
                    '%s: Value for option %s must be a single logical value.\n',...
                    app,...
                    name);
                return;
            end
            saveAllFields = value;
        case 'outfile'
            if ~ischar(value)
                fprintf(2,...
                    '%s: Value for option %s must be a string indicating a file name for output.\n',...
                    app,...
                    name);
                return;
            end
            if ~strcmp(value(end-2:end),'.nc')
                value=[value '.nc'];
            end
            outFile = value;
        case 'maxseparation'
            if numel(value)~=1 | ~isnumeric(value)
                fprintf(2,...
                    '%s: Value for option %s must be a single numeric value in kilometers.\n',...
                    app,...
                    name);
                return;
            end
            maxSeparation = value;
        case 'totalscomparisonmethod'
            if ~ischar(value) | ~ismember(lower(value),totalsMethodOptions)
                methodlist=char(strcat(totalsMethodOptions,','))';
                methodlist=methodlist(1:end-1);
                fprintf(2,...
                    '%s: Value for option %s must be a one of the following: %s\n',...
                    app,...
                    name,...
                    methodlist);
                return;
            end
            totalsComparisonMethod = value;
        case 'totalstoremove'
            if ~ischar(value)
                fprintf(2,...
                    '%s: Value for option %s must be a string indicating expression to use for bad totals to remove before matching.\n',...
                    app,...
                    name);
                return;
            end
            totalsToRemove = value;
        case 'radialstoremove'
            if ~ischar(value)
                fprintf(2,...
                    '%s: Value for option %s must be a string indicating expression to use for bad radials to remove before matching.\n',...
                    app,...
                    name);
                return;
            end
            radialsToRemove = value;
        otherwise
            fprintf(2,...
                '%s: Invalid option specified: %s.\n',...
                app,...
                name);
    end
end


if ischar(processedDrifterData) & exist(processedDrifterData,'file')==2
    time=datenum(1970,1,1,0,0,double(ncread(processedDrifterData,'time')));
    lon=ncread(processedDrifterData,'lon');
    lat=ncread(processedDrifterData,'lat');

    drifterVelocities.time=time;
    drifterVelocities.lon=lon;
    drifterVelocities.lat=lat;
    drifterVelocities.u=ncread(processedDrifterData,'u')*100;
    drifterVelocities.v=ncread(processedDrifterData,'v')*100;
    drifterVelocities.units='cm/s';
elseif isstruct(processedDrifterData)
    drifterVelocities=processedDrifterData;
    time=drifterVelocities.time;
    lon=drifterVelocities.lon;
    lat=drifterVelocities.lat;
    drifterVelocities.u=drifterVelocities.u*100;
    drifterVelocities.v=drifterVelocities.v*100;
    drifterVelocities.units='cm/s';
else
    fprintf(2,...
        '%s: first input to function must be a netcdf file containing processed drifter data or a structured array as output from processDrifterFiles.m.\n',...
        app);
    return;
end

radialSpeeds=struct([]);
totalsVelocities=struct([]);


if ~strcmp(lower(radialType),'none')
    % all site options
    allRadials=dir(radialDir);
    allRadials=struct2cell(allRadials);
    ind=find(~strncmp(allRadials(1,:),'.',1));
    allRadials=allRadials(1,ind);
    if numel(radialSites)==1 & strncmpi(radialSites{1},'all',3)
        if exist('HFR_sites.xlsx')
            [~,~,siteInfo]=xlsread('HFR_sites.xlsx');
            lonC=find(strcmpi('lon',siteInfo(1,:)));
            latC=find(strcmpi('lat',siteInfo(1,:)));
            rangeC=find(strcmpi('range',siteInfo(1,:)));
            siteC=find(strcmpi('site',siteInfo(1,:)));
            siteTypes=radialSites{1};
            radialSites=[];
            if(ismember(lower(siteTypes),{'all','all5','all5mhz','alllong','alllr'}))
                ind=find(strcmpi(siteInfo(:,rangeC),'long'));
                for n=1:length(ind)
                    DI = distance(lat,lon,...
                        siteInfo{ind(n),latC}*ones(size(lat)),siteInfo{ind(n),lonC}*ones(size(lon)));
                    DI = deg2km(DI);
                    if sum(DI<100)>min(25,length(lat)/2)
                        radialSites=[radialSites,siteInfo(ind(n),siteC)];
                    end
                end
            end
            if(ismember(lower(siteTypes),{'all','all13','all13mhz','allmedium','allmr'}))
                ind=find(strcmpi(siteInfo(:,rangeC),'medium'));
                for n=1:length(ind)
                    DI = distance(lat,lon,...
                        siteInfo{ind(n),latC}*ones(size(lat)),siteInfo{ind(n),lonC}*ones(size(lon)));
                    DI = deg2km(DI);
                    if sum(DI<75)>min(25,length(lat)/2)
                        radialSites=[radialSites,siteInfo(ind(n),siteC)];
                    end
                end
            end
            if(ismember(lower(siteTypes),{'all','all25','all25mhz','allstandard','allshort','allsr'}))
                ind=find(strcmpi(siteInfo(:,rangeC),'standard'));
                for n=1:length(ind)
                    DI = distance(lat,lon,...
                        siteInfo{ind(n),latC}*ones(size(lat)),siteInfo{ind(n),lonC}*ones(size(lon)));
                    DI = deg2km(DI);
                    if sum(DI<50)>min(25,length(lat)/2)
                        radialSites=[radialSites,siteInfo(ind(n),siteC)];
                    end
                end
            end
        else
            radialSites=allRadials;
        end
    end
    
    missingSites=setdiff(radialSites,allRadials);
    radialSites=intersect(radialSites,allRadials);
    
    if ~isempty(missingSites)
        missingSites=char(strcat(missingSites,','))';
        missingSites=missingSites(1:end-1);
        warning('%s: The following radial sites are not available in directory %s and will be skipped: %s\n',...
            app,...
            radialDir,...
            missingSites);
    end
    
    for s=1:length(radialSites)
        site=radialSites{s};
        radialSpeeds(1).attributes.radial_directory=radialDir;
        radialSpeeds.attributes.furthest_valid_distance=maxSeparation;
        radialSpeeds.attributes.furthest_valid_distance_units='km';
        radialSpeeds.attributes.radials_excluded=radialsToRemove{1};
        for ra=2:length(radialsToRemove)
            radialSpeeds.attributes.radials_excluded=...
                [radialSpeeds.attributes.radials_excluded, 'or', ...
                radialsToRemove{ra}];
        end
        if ismember(lower(radialType),{'both','measured','rdlm'})
            type='measured';
            radialSpeeds.(site).(type).time=time;
            for ti=1:length(time)
                t=time(ti);
                if isdir([radialDir site datestr(t,'/yyyy_mm/')])
                    radialFile=dir([radialDir, site, datestr(t,'/yyyy_mm/'),...
                        'RDLm_' site '_' datestr(t,'yyyy_mm_dd_HH00') '*']);
                    radialDirTemp=[radialDir, site, datestr(t,'/yyyy_mm/')];
                else
                    radialFile=dir([radialDir, site, ...
                        'RDLm_' site '_' datestr(t,'yyyy_mm_dd_HH00') '*']);
                    radialDirTemp=[radialDir, site];
                end
                if ~isempty(radialFile)
                    radialFile=[radialDirTemp radialFile.name];
                    [closestRadial,siteOrigin]=getNearestRadial(radialFile,...
                        drifterVelocities.lon(ti),...
                        drifterVelocities.lat(ti),...
                        drifterVelocities.u(ti),...
                        drifterVelocities.v(ti),...
                        maxSeparation,...
                        radialsToRemove);
                    newFields=fields(closestRadial);
                    if(~isfield(radialSpeeds.(site).(type),'SiteOrigin')|...
                        isempty(radialSpeeds.(site).(type).SiteOrigin))
                        radialSpeeds.(site).(type).SiteOrigin=siteOrigin;
                    end
                    for f=1:length(newFields)
                        if ~isfield(radialSpeeds.(site).(type),newFields{f})
                            radialSpeeds.(site).(type).(newFields{f})=...
                                nan(size(radialSpeeds.(site).(type).time));
                        end
                        radialSpeeds.(site).(type).(newFields{f})(ti)=closestRadial.(newFields{f});
                    end
                else
                    warning('%s: No %s radials found for %s at time %s\n',...
                        app,...
                        type,...
                        site,...
                        datestr(t,'yyyy-mm-dd HH:00'));
                end
            end
        end
        
        if ismember(lower(radialType),{'both','ideal','rdli'})
            type='ideal';
            radialSpeeds.(site).(type).time=time;
            for ti=1:length(time)
                t=time(ti);
                if isdir([radialDir site datestr(t,'/yyyy_mm/')])
                    radialFile=dir([radialDir, site, datestr(t,'/yyyy_mm/'),...
                        'RDLi_' site '_' datestr(t,'yyyy_mm_dd_HH00') '*']);
                    radialDirTemp=[radialDir, site, datestr(t,'/yyyy_mm/')];
                else
                    radialFile=dir([radialDir, site, ...
                        'RDLi_' site '_' datestr(t,'yyyy_mm_dd_HH00') '*']);
                    radialDirTemp=[radialDir, site];
                end
                if ~isempty(radialFile)
                    radialFile=[radialDirTemp radialFile.name];
                    closestRadial=getNearestRadial(radialFile,...
                        drifterVelocities.lon(ti),...
                        drifterVelocities.lat(ti),...
                        drifterVelocities.u(ti),...
                        drifterVelocities.v(ti),...
                        maxSeparation,...
                        radialsToRemove);
                    newFields=fields(closestRadial);
                    for f=1:length(newFields)
                        if ~isfield(radialSpeeds.(site).(type),newFields{f})
                            radialSpeeds.(site).(type).(newFields{f})=...
                                nan(size(radialSpeeds.(site).(type).time));
                        end
                        radialSpeeds.(site).(type).(newFields{f})(ti)=closestRadial.(newFields{f});
                    end
                else
                    warning('%s: No %s radials found for %s at time %s\n',...
                        app,...
                        type,...
                        site,...
                        datestr(t,'yyyy-mm-dd HH:00'));
                end
            end
        end
    end
end

if ~strcmp(totalsNetworks{1},'none')
    if strcmp(totalsComparisonMethod,'nearest')
        sx=maxSeparation;
        sy=maxSeparation;
        type='nearest';
    elseif strcmp(totalsComparisonMethod,'mean')
        sx=maxSeparation;
        sy=maxSeparation;
        type='mean';
    elseif strcmp(totalsComparisonMethod,'median')
        sx=maxSeparation;
        sy=maxSeparation;
        type='median';
    elseif ismember(totalsComparisonMethod,{'weightedgaussian','weightedgaussianellipse'})
        sx=maxSeparation;
        sy=maxSeparation;
        type='gaussian';
    elseif ismember(totalsComparisonMethod,{'weightedexponential','weightedexponentialellipse'})
        sx=maxSeparation;
        sy=maxSeparation;
        type='exponential';
    end
    if ~strcmp(totalsDir,'thredds')
        totalsVelocities(1).attributes.totals_directory=totalsDir;
        totalsVelocities.attributes.furthest_valid_distance=maxSeparation;
        totalsVelocities.attributes.furthest_valid_distance_units='km';
        totalsVelocities.attributes.totals_excluded=totalsToRemove{1};
        for ta=2:length(totalsToRemove)
            totalsVelocities.attributes.totals_excluded=...
                [totalsVelocities.attributes.totals_excluded, 'or', ...
                totalsToRemove{ta}];
        end
        totalsVelocities.attributes.totals_matching_method=type;
        if ismember(totalsComparisonMethod,{'weightedexponential','weightedgaussian'})
            totalsVelocities.attributes.totals_matching_method=['weighted ' type];
        elseif ismember(totalsComparisonMethod,{'weightedexponentialellipse','weightedgaussianellipse'})
            totalsVelocities.attributes.totals_matching_method=['weighted ' type ' ellipse'];
        end
        totalsVelocities.totals.time=time;
        for ti=1:length(time)
            t=time(ti);
            disp(['totals ' datestr(t)])
            if isdir([totalsDir datestr(t,'yyyy_mm/')])
                totalsFile=dir([totalsDir, datestr(t,'yyyy_mm/'),...
                    '*' datestr(t,'yyyy_mm_dd_HH00') '*.nc']);
                totalsDirTemp=[totalsDir, datestr(t,'yyyy_mm/')];
            else
                totalsFile=dir([totalsDir, ...
                    '*' datestr(t,'yyyy_mm_dd_HH00') '*.nc']);
                totalsDirTemp=totalsDir;
            end
            if isempty(totalsFile)
                warning('%s: No totals found at time %s in provided directory\n',...
                    app,...
                    datestr(t,'yyyy-mm-dd HH:00'));
            elseif length(totalsFile)>1
                warning('%s: Multiple totals found at time %s in provided directory. Skipping.\n',...
                    app,...
                    datestr(t,'yyyy-mm-dd HH:00'));
            elseif totalsFile.bytes<80000
                warning('%s: File %s is very small and may be corrupt. Skipping.\n',...
                    app,...
                    totalsFile.name);
            else
                totalsFile=[totalsDirTemp totalsFile.name];
                if ismember(totalsComparisonMethod,{'weightedexponentialellipse','weightedgaussianellipse'})
                    i1=max(ti-1,1);
                    i2=min(ti+1,length(time));
                    DI = distance(lat(i1:i2-1),lon(i1:i2-1),lat(i1+1:i2),lon(i1+1:i2));
                    DI = deg2km(DI);
                    AZ =  azimuth(lat(i1:i2-1),lon(i1:i2-1),lat(i1+1:i2),lon(i1+1:i2));
                    AZ = (90-AZ)*pi/180;
                    s=[nanmean(cos(AZ).*DI),nanmean(sin(AZ).*DI)];
                    sx=s(1)*maxSeparation/max(s);
                    sy=s(2)*maxSeparation/max(s);
                end
                try
                    totalsData.lon=ncread(totalsFile,'lon');
                    totalsData.lat=ncread(totalsFile,'lat');
                    totalsData.u=ncread(totalsFile,'u');
                    totalsData.v=ncread(totalsFile,'v');
                    totalsData.u_err=ncread(totalsFile,'u_err');
                    totalsData.v_err=ncread(totalsFile,'v_err');
                    totalsData.num_radials=ncread(totalsFile,'num_radials');
                    matchingTotal=getMatchingTotal(totalsData,...
                        drifterVelocities.lon(ti),...
                        drifterVelocities.lat(ti),...
                        sx,sy,type,...
                        totalsToRemove);
                    newFields=fields(matchingTotal);
                    for f=1:length(newFields)
                        if ~isfield(totalsVelocities.totals,newFields{f})
                            totalsVelocities.totals.(newFields{f})=...
                                nan(size(totalsVelocities.totals.time));
                        end
                        totalsVelocities.totals.(newFields{f})(ti)=matchingTotal.(newFields{f});
                    end
                catch
                    warning('%s: Issue reading from file %s\n',...
                        app,...
                        totalsFile);
                end
            end
        end
    else
        totalsVelocities(1).attributes.totals_directory='http://tds.marine.rutgers.edu/thredds/cool/codar/cat_totals.html';
        totalsVelocities.attributes.furthest_valid_distance=maxSeparation;
        totalsVelocities.attributes.furthest_valid_distance_units='km';
        totalsVelocities.attributes.totals_excluded=totalsToRemove{1};
        for ta=2:length(totalsToRemove)
            totalsVelocities.attributes.totals_excluded=...
                [totalsVelocities.attributes.totals_excluded, 'or', ...
                totalsToRemove{ta}];
        end
        totalsVelocities.attributes.totals_matching_method=type;
        if ismember(totalsComparisonMethod,{'weightedexponential','weightedgaussian'})
            totalsVelocities.attributes.totals_matching_method=['weighted ' type];
        elseif ismember(totalsComparisonMethod,{'weightedexponentialellipse','weightedgaussianellipse'})
            totalsVelocities.attributes.totals_matching_method=['weighted ' type ' ellipse'];
        end
        for network=1:length(totalsNetworks)
            if ismember(lower(totalsNetworks{network}),...
                    {'all','25','25mhz','sr','standardrange','shortrange'})
                totalsFile='http://tds.marine.rutgers.edu/thredds/dodsC/cool/codar/totals/25Mhz_1km_realtime_fmrc/Maracoos_25MHz_1km_Totals-FMRC_best.ncd';
                totalsData.lon=ncread(totalsFile,'lon');
                totalsData.lat=ncread(totalsFile,'lat');
                if sum(drifterVelocities.lon>=min(totalsData.lon) & ...
                        drifterVelocities.lon<=max(totalsData.lon) & ...
                        drifterVelocities.lat>=min(totalsData.lat) & ...
                        drifterVelocities.lat<=max(totalsData.lat)) > ...
                        min(50,length(time)/10)
                    totalsVelocities.totals_25MHz.time=time;
                    totalsTime=datenum(2012,8,25,double(ncread(totalsFile,'time'))+17,0,0);
                    for ti=1:length(time)
                        t=time(ti);
                        t_ind=find(abs(totalsTime-t)<.01);
                        if length(t_ind)==1
                            if ismember(totalsComparisonMethod,{'weightedexponentialellipse','weightedgaussianellipse'})
                                i1=max(ti-1,1);
                                i2=min(ti+1,length(time));
                                DI = distance(lat(i1:i2-1),lon(i1:i2-1),lat(i1+1:i2),lon(i1+1:i2));
                                DI = deg2km(DI);
                                AZ =  azimuth(lat(i1:i2-1),lon(i1:i2-1),lat(i1+1:i2),lon(i1+1:i2));
                                AZ = (90-AZ)*pi/180;
                                s=[nanmean(cos(AZ)*DI),nanmean(sin(AZ)*DI)];
                                sx=s(1)*maxSeparation/max(s);
                                sy=s(2)*maxSeparation/max(s);
                            end
                            totalsData.u=ncread(totalsFile,'u',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.v=ncread(totalsFile,'v',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.u_err=ncread(totalsFile,'u_err',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.v_err=ncread(totalsFile,'v_err',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.num_radials=ncread(totalsFile,'num_radials',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            matchingTotal=getMatchingTotal(totalsData,...
                                drifterVelocities.lon(ti),...
                                drifterVelocities.lat(ti),...
                                sx,sy,type,...
                                totalsToRemove);
                            newFields=fields(matchingTotal);
                            for f=1:length(newFields)
                                if ~isfield(totalsVelocities.totals_25MHz,newFields{f})
                                    totalsVelocities.totals_25MHz.(newFields{f})=...
                                        nan(size(totalsVelocities.totals_25MHz.time));
                                end
                                totalsVelocities.totals_25MHz.(newFields{f})(ti)=matchingTotal.(newFields{f});
                            end
                        end
                    end
                end
            end
            
            if ismember(lower(totalsNetworks{network}),...
                    {'all','13','13mhz','mr','mediumrange'})
                totalsFile='http://tds.marine.rutgers.edu/thredds/dodsC/cool/codar/totals/13Mhz_2km_realtime_fmrc/Maracoos_13MHz_2km_Totals-FMRC_best.ncd';
                totalsData.lon=ncread(totalsFile,'lon');
                totalsData.lat=ncread(totalsFile,'lat');
                if sum(drifterVelocities.lon>=min(totalsData.lon) & ...
                        drifterVelocities.lon<=max(totalsData.lon) & ...
                        drifterVelocities.lat>=min(totalsData.lat) & ...
                        drifterVelocities.lat<=max(totalsData.lat)) > ...
                        min(50,length(time)/10)
                    totalsVelocities.totals_13MHz.time=time;
                    totalsTime=datenum(2011,8,1,double(ncread(totalsFile,'time')),0,0);
                    for ti=1:length(time)
                        t=time(ti);
                        t_ind=find(abs(totalsTime-t)<.01);
                        if length(t_ind)==1
                            if ismember(totalsComparisonMethod,{'weightedexponentialellipse','weightedgaussianellipse'})
                                i1=max(ti-1,1);
                                i2=min(ti+1,length(time));
                                DI = distance(lat(i1:i2-1),lon(i1:i2-1),lat(i1+1:i2),lon(i1+1:i2));
                                DI = deg2km(DI);
                                AZ =  azimuth(lat(i1:i2-1),lon(i1:i2-1),lat(i1+1:i2),lon(i1+1:i2));
                                AZ = (90-AZ)*pi/180;
                                s=[nanmean(cos(AZ)*DI),nanmean(sin(AZ)*DI)];
                                sx=s(1)*maxSeparation/max(s);
                                sy=s(2)*maxSeparation/max(s);
                            end
                            totalsData.u=ncread(totalsFile,'u',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.v=ncread(totalsFile,'v',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.u_err=ncread(totalsFile,'u_err',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.v_err=ncread(totalsFile,'v_err',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.num_radials=ncread(totalsFile,'num_radials',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            matchingTotal=getMatchingTotal(totalsData,...
                                drifterVelocities.lon(ti),...
                                drifterVelocities.lat(ti),...
                                sx,sy,type,...
                                totalsToRemove);
                            newFields=fields(matchingTotal);
                            for f=1:length(newFields)
                                if ~isfield(totalsVelocities.totals_13MHz,newFields{f})
                                    totalsVelocities.totals_13MHz.(newFields{f})=...
                                        nan(size(totalsVelocities.totals_13MHz.time));
                                end
                                totalsVelocities.totals_13MHz.(newFields{f})(ti)=matchingTotal.(newFields{f});
                            end
                        end
                    end
                end
            end
            
            if ismember(lower(totalsNetworks{network}),...
                    {'all','5','5mhz','lr','longrange'})
                totalsFile='http://tds.marine.rutgers.edu/thredds/dodsC/cool/codar/totals/5Mhz_6km_realtime_fmrc/Maracoos_5MHz_6km_Totals-FMRC_best.ncd';
                totalsData.lon=ncread(totalsFile,'lon');
                totalsData.lat=ncread(totalsFile,'lat');
                if sum(drifterVelocities.lon>=min(totalsData.lon) & ...
                        drifterVelocities.lon<=max(totalsData.lon) & ...
                        drifterVelocities.lat>=min(totalsData.lat) & ...
                        drifterVelocities.lat<=max(totalsData.lat)) > ...
                        min(50,length(time)/10)
                    totalsVelocities.totals_5MHz.time=time;
                    totalsTime=datenum(2006,1,1,double(ncread(totalsFile,'time')),0,0);
                    for ti=1:length(time)
                        t=time(ti);
                        t_ind=find(abs(totalsTime-t)<.01);
                        if length(t_ind)==1
                            if ismember(totalsComparisonMethod,{'weightedexponentialellipse','weightedgaussianellipse'})
                                i1=max(ti-1,1);
                                i2=min(ti+1,length(time));
                                DI = distance(lat(i1:i2-1),lon(i1:i2-1),lat(i1+1:i2),lon(i1+1:i2));
                                DI = deg2km(DI);
                                AZ =  azimuth(lat(i1:i2-1),lon(i1:i2-1),lat(i1+1:i2),lon(i1+1:i2));
                                AZ = (90-AZ)*pi/180;
                                s=[nanmean(cos(AZ)*DI),nanmean(sin(AZ)*DI)];
                                sx=s(1)*maxSeparation/max(s);
                                sy=s(2)*maxSeparation/max(s);
                            end
                            totalsData.u=ncread(totalsFile,'u',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.v=ncread(totalsFile,'v',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.u_err=ncread(totalsFile,'u_err',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.v_err=ncread(totalsFile,'v_err',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            totalsData.num_radials=ncread(totalsFile,'num_radials',[1 1 t_ind],[length(totalsData.lon) length(totalsData.lat) 1]);
                            matchingTotal=getMatchingTotal(totalsData,...
                                drifterVelocities.lon(ti),...
                                drifterVelocities.lat(ti),...
                                sx,sy,type,...
                                totalsToRemove);
                            newFields=fields(matchingTotal);
                            for f=1:length(newFields)
                                if ~isfield(totalsVelocities.totals_5MHz,newFields{f})
                                    totalsVelocities.totals_5MHz.(newFields{f})=...
                                        nan(size(totalsVelocities.totals_5MHz.time));
                                end
                                totalsVelocities.totals_5MHz.(newFields{f})(ti)=matchingTotal.(newFields{f});
                            end
                        end
                    end
                end
            end
        end
    end
end


if ~strcmp(outFile,'none')
    drifterHFRmatches2nc(processedDrifterData, outFile, saveAllFields, ...
        drifterVelocities, radialSpeeds, totalsVelocities);
end