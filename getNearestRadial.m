function [closestRadial,siteOrigin]=getNearestRadial(radialFile,...
    drifterlon, drifterlat,...
    drifteru, drifterv,...
    maxSeparation, radialsToRemove)

closestRadial.rotated_drifter_velocity=nan;
closestRadial.distance_to_HFR_radial=nan;
closestRadial.HFR_radial_velocity=nan;
closestRadial.HFR_radial_lon=nan;
closestRadial.HFR_radial_lat=nan;

rdl=loadRDLFile(radialFile,true);
siteOrigin=rdl.SiteOrigin;

if isempty(rdl.LonLat)
    warning(['No data in ' radialFile])
    return;
end

% find index of table-start row in header
tablestart_ind=find(strncmp('%TableStart',rdl.OtherMetadata.Header,length('%TableStart')));
tablestart_ind=tablestart_ind(1);

% divide into pre-data and post-data subsets of header
subset1=rdl.OtherMetadata.Header(1:tablestart_ind+2);
%subset3=rdl.OtherMetadata.Header(tablestart_ind+3:end);

% get number of columns (in header description)
colnum_ind=find(strncmp('%TableColumns',subset1,length('%TableColumns')));
colnum=str2num(subset1{colnum_ind}(length('%TableColumns:  '):end));

% get 4-char variable IDs
coltype_ind=strncmp('%TableColumnTypes',subset1,length('%TableColumnTypes'));
coltype=subset1{coltype_ind};
ind=find(coltype==':');
coltype=coltype(ind+1:end);
coltype(coltype==' ')='';
if(length(coltype)~=colnum*4)
    warning('column type labels do not match up with column count')
end

% get long-name header to each variable
header=cell(1,colnum);
for k=1:colnum
    header{k}=coltype(k*4-3:k*4);
end

for k=1:length(header)
    eval([header{k} '=rdl.OtherMetadata.RawData(:,k);'])
end

% remove any bad data as indicated in radialsToRemove
ind_bad=[];
if ~strcmp(radialsToRemove{1},'none')
    for rqc=1:length(radialsToRemove)
        try
            eval(['ind_bad_new=find(' radialsToRemove{rqc} ');'])
            ind_bad=union(ind_bad,ind_bad_new);
        end
    end
end

DI = distance(LATD,LOND,drifterlat*ones(size(LATD)),drifterlon*ones(size(LOND)));
DI = deg2km(DI);
ind_bad_new=find(DI>maxSeparation);
ind_bad=union(ind_bad,ind_bad_new);

DI(ind_bad)=[];
for k=1:length(header)
    eval([header{k} '(ind_bad)=[];'])
end

if ~isempty(DI)
    [~,ind_closest]=min(DI);
    add_fields=setdiff(header,{'LOND','LATD','VELU','VELV','VELO'});
    closestRadial.distance_to_HFR_radial=DI(ind_closest);
    closestRadial.HFR_radial_velocity=VELO(ind_closest);
    closestRadial.HFR_radial_lon=LOND(ind_closest);
    closestRadial.HFR_radial_lat=LATD(ind_closest);
    for k=1:length(add_fields)
        eval(['closestRadial.HFR_radial_', add_fields{k}, ...
            '=', add_fields{k} '(ind_closest);'])
    end
    closestRadial.rotated_drifter_velocity=-rot(drifteru,drifterv,...
        closestRadial.HFR_radial_BEAR);
end

    