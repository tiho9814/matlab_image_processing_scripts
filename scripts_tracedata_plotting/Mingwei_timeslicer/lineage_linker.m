function raw = lineage_linker(cellid, cells)
% this function finds a cells progenitors and progenies and link all their
% tracedata togethter. For progenies, it randomly choose one daughter
% between the two.
% use cellinfo input.
% Mingwei Min 2017.4.22

raw         = cells(cellid).raw;

% find ancestors 
condition   = true;
tempid      = cellid;
while condition
    if ~(cells(tempid).mother==0)
        tempid      = cells(tempid).mother;
        goodframes  = find(~isnan(cells(tempid).raw(:,1)));
        if ~cells(tempid).noisy
            raw(goodframes,:)   = cells(tempid).raw(goodframes,:);
        end
    else
        condition = false;
    end
end

% find progenies (randomly pick a daughter between the two)
condition   = true;
tempid      = cellid;
while condition
    if ~(cells(tempid).daughter==0)
        tempid      = cells(tempid).daughter(1);
        goodframes  = find(~isnan(cells(tempid).raw(:,1)));
        if ~cells(tempid).noisy
            raw(goodframes,:)   = cells(tempid).raw(goodframes,:);
        end
    else
        condition = false;
    end
end