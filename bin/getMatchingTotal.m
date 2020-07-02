function matchingTotal=getMatchingTotal(totalsData,...
    drifterlon,drifterlat,...
    sx,sy,type,...
    totalsToRemove)


matchingTotal.HFR_totals_u=nan;
matchingTotal.HFR_totals_v=nan;

all_fields=fields(totalsData);
for f=1:length(all_fields)
    eval([all_fields{f} '=totalsData.' all_fields{f} ';'])
end

ind_bad=[];
if ~strcmp(totalsToRemove{1},'none')
    for rqc=1:length(totalsToRemove)
        try
            eval(['ind_bad_new=find(' totalsToRemove{rqc} ');'])
            ind_bad=union(ind_bad,ind_bad_new);
        end
    end
end

[LATD,LOND]=meshgrid(lat,lon);
DI = distance(LATD,LOND,drifterlat*ones(size(LATD)),drifterlon*ones(size(LOND)));
DI = deg2km(DI);
ind_bad_new=find(DI>max(sx,sy));
ind_bad=union(ind_bad,ind_bad_new);

u(ind_bad)=nan;
v(ind_bad)=nan;
DI(ind_bad)=nan;
DI(isnan(u)|isnan(v))=nan;

if ~all(isnan(DI(:)))
    if strcmp(type,'nearest')
        [~,ind_closest]=min(DI(:));
        matchingTotal.HFR_totals_u=u(ind_closest);
        matchingTotal.HFR_totals_v=v(ind_closest);
        matchingTotal.distance_to_closest_HFR_total=DI(ind_closest);
        matchingTotal.HFR_totals_u_err=u_err(ind_closest);
        matchingTotal.HFR_totals_v_err=v_err(ind_closest);
        matchingTotal.HFR_totals_num_radials=num_radials(ind_closest);
    elseif ismember(type,{'mean','median'})
        eval(['matchingTotal.HFR_totals_u=nan' type '(u(:));'])
        eval(['matchingTotal.HFR_totals_v=nan' type '(v(:));'])
        matchingTotal.HFR_num_totals=sum(~isnan(DI(:)));
    else
        AZ = azimuth(LATD,LOND,drifterlat*ones(size(LATD)),drifterlon*ones(size(LOND)));
        AZ = (90-AZ)*pi/180;
        dx=cos(AZ).*DI;
        dy=sin(AZ).*DI;
        if strcmp(type,'gaussian')
            wt=exp(-(dx.^2/sx^2 + dy.^2/sy^2));
        elseif strcmp(type,'exponential')
            wt=exp(-sqrt(dx.^2/sx^2 + dy.^2/sy^2));
        end
        u(isnan(DI))=nan;
        v(isnan(DI))=nan;
        wt(isnan(DI))=nan;
        matchingTotal.HFR_totals_u=nansum(u(:).*wt(:))/nansum(wt(:));
        matchingTotal.HFR_totals_v=nansum(v(:).*wt(:))/nansum(wt(:));
        matchingTotal.HFR_num_totals=sum(~isnan(DI(:)));
    end
end

