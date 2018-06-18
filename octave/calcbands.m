function eng = calcbands(kptlist)
%
% This function calculates the energies on a list of kpoints.
%
% global variables:
%   orb: fractional coordinates of local orbitals
%
% input:
%   kptlist: list of fractional coordinates of kpoints, nkpt*3 matrix
%
% output:
%   eng: calcualted energies, nkpt*norb matrix
%
    global orb;
    
    [nkpt, foo] = size(kptlist);
    [norb, foo] = size(orb);
    
    eng = zeros(nkpt, norb);
    for ikpt = 1:nkpt
        Hk = setham(kptlist(ikpt,:));
        eng(ikpt,:) = sort(real(eig(Hk)))';
    end
