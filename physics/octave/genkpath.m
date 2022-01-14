function kpath = genkpath(kptlist, ninterp)
%
% This function generates the kpath that connects the given highly symmetric
% kpoints in the first Brillouin zone.
%
% input:
%   kptlist: coordinates of highly symmetric k-points, N*3 array, 
%            does not need to be fractional
%   ninterp: numbers of interpolated k-points on each segment of k-path, integer
%
% output:
%   kpath: coordinates of generated k-path, array with 3 columns
%
    [nkpt, foo] = size(kptlist);
    kpath = zeros((nkpt-1)*ninterp+1, 3);
    icount = 1;
    for ikpt = 1:nkpt-1
        k0 = kptlist(ikpt,:);
        k1 = kptlist(ikpt+1,:);
        for iint = 1:ninterp
            kpath(icount,:) = k0 + (k1-k0) * (iint-1) / ninterp;
            icount = icount + 1;
        end
    end
    kpath(icount,:) = kptlist(nkpt,:);
