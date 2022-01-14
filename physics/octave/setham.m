function ham = setham(kpt)
%
% This function builds the Hamiltonian at a given kpoint.
%
% global variables:
%   orb: fractional coordinates of local orbitals
%   ons: on-site energies, with norb elements, either real or complex, should be
%        sorted according to local orbital indices
%   hop: hopping integrals, with 6 columns, 1-2 for i/j, 3-5 for Rn and 6 for tij,
%        tij may be either real or complex
%   set_conj: whether to set conjugate terms automatically
%
% input:
%   kpt: fractional coordinates of the single k-point, 1*3 array
%
% output:
%   ham: the Hamiltonian, square matrix of norb*norb
%
% This is a core-level function called by other wrappers.
% In most cases you do not need to access it directly.
%
    global orb ons hop set_conj;

    [norb, foo] = size(orb);
    ham = zeros(norb, norb);

    for iorb = 1:norb
        ham(iorb,iorb) = ons(iorb);
    end

    [nhop, foo] = size(hop);
    if (set_conj == true)
        for ihop = 1:nhop
            orbi = int32(hop(ihop,1));
            orbj = int32(hop(ihop,2));
            Rij = orb(orbj,:) + hop(ihop,3:5) - orb(orbi,:);
            tij = hop(ihop,6);
            hij = tij * exp(i * 2*pi * kpt * Rij');
            ham(orbi,orbj) = ham(orbi,orbj) + hij;
            ham(orbj,orbi) = ham(orbj,orbi) + conj(hij);
        end
    else
        for ihop = 1:nhop
            orbi = int32(hop(ihop,1));
            orbj = int32(hop(ihop,2));
            Rij = orb(orbj,:) + hop(ihop,3:5) - orb(orbi,:);
            tij = hop(ihop,6);
            hij = tij * exp(i * 2*pi * kpt * Rij');
            ham(orbi,orbj) = ham(orbi,orbj) + hij;
        end
    end
