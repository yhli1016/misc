function [eng, proj] = projbands(kptlist, lolist)
%
% This function calculates the energies and projections
% on given local orbitals on a list of kpoints.
%
% global variables:
%   orb: fractional coordinates of local orbitals
%
% input:
%   kptlist: list of fractional coordinates of kpoints, nkpt*3 matrix
%   lolist: list of indices of orbitals to project on
%
% output:
%   eng: calcualted energies, nkpt*norb matrix
%   proj: projection on given local orbitals, nkpt*borb matrix
%
    global orb;
    
    [nkpt, foo] = size(kptlist);
    [norb, foo] = size(orb);
    nlo = length(lolist);
    
    eng = zeros(nkpt, norb);
    proj = zeros(nkpt, norb);

    % loop over all kpoints
    for ikpt = 1:nkpt
        % build H(k) and diagonalize it
        Hk = setham(kptlist(ikpt,:));
        [eigvec, eigval] = eig(Hk);
        
        % loop over all eigvenstates to collect energies and projections
        for iorb = 1:norb
            % eigval is a diagonal matrix
            eng(ikpt,iorb) = real(eigval(iorb,iorb));

            % each column of eigvec is the corresponding eigenvector
            for ilo = 1:nlo
                proj(ikpt,iorb) = proj(ikpt,iorb) + abs(eigvec(lolist(ilo),iorb)) ^ 2;
            end
            % normalize
            proj(ikpt,iorb) = proj(ikpt,iorb) / sum(abs(eigvec(:,iorb)) .^ 2);
        end
    end
