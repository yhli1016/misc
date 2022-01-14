function bm(dat, nband_val)
% This function calculates VBM / CBM / Indirect Eg / Direct Eg from band data.
%
% dat: band data which must be in the following form
%   kpath(1) band1(1) band2(1) band3(1) ... bandm(1)
%   kpath(2) band1(2) band2(2) band3(2) ... bandm(2)
%   ......            ......
%   kpath(n) band1(n) band2(n) band3(n) ... bandm(n)
%
% nband_val: number of valence bands in band data

% split band data
[nk, ncol] = size(dat);
band_full = dat(:, 2:ncol);
[nk, nb] = size(band_full);
band_val = band_full(:, 1:nband_val);
band_con = band_full(:, nband_val+1:nb);

% determine vbmax and cbmin for each kpoint
vbmax = zeros(nk, 2);
cbmin = zeros(nk, 2);
for ik = 1:nk
    [eng, ib] = max(band_val(ik,:));
    vbmax(ik,:) = [eng, ib];
    [eng, ib] = min(band_con(ik,:));
    cbmin(ik,:) = [eng, ib];
end

% determine indirect band gap
[VBM, ik_VBM] = max(vbmax(:,1));
[CBM, ik_CBM] = min(cbmin(:,1));
IEg = CBM - VBM;
ib_VBM = vbmax(ik_VBM, 2);
ib_CBM = cbmin(ik_CBM, 2);
printf('VBM = %10.5f\n', VBM);
printf('CBM = %10.5f\n', CBM);
printf('IEg = %10.5f, (%4d,%4d) -> (%4d,%4d)\n', IEg, ib_VBM, ik_VBM, ib_CBM, ik_CBM);

% determine direct band gap
deltaE = cbmin(:,1) - vbmax(:,1);
[DEg, ik_DEg] = min(deltaE);
ib_DEgV = vbmax(ik_DEg, 2);
ib_DEgC = cbmin(ik_DEg, 2);
printf('DEg = %10.5f, (%4d,%4d) -> (%4d,%4d)\n', DEg, ib_DEgV, ik_DEg, ib_DEgC, ik_DEg);
