## Copyright (C) 2017 yhli
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} uniformk (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: yhli <yhli@PC05>
## Created: 2017-05-01

function uniformk (meshsize)
%
% This function generates uniform kgrid in the first Brillouin zone and write it
% to 'kgrid.dat'.
%
% meshsize: size of the kgrid, a vector containing 3 integers
% kgrid: thus generated list of kpoints, with meshsize(1) * meshsize(2) *
%        meshsize(3) rows and 3 columns
%
    nk1 = meshsize(1);
    nk2 = meshsize(2);
    nk3 = meshsize(3);
    xv1 = linspace(0, (nk1-1)/nk1, nk1);
    xv2 = linspace(0, (nk2-1)/nk2, nk2);
    xv3 = linspace(0, (nk3-1)/nk3, nk3);
    kgrid = zeros(nk1*nk2*nk3, 3);
    icount = 1;
    for ik1 = 1:nk1
        for ik2 = 1:nk2
            for ik3 = 1:nk3
                kgrid(icount,:) = [xv1(ik1), xv2(ik2), xv3(ik3)];
                icount = icount + 1;
            endfor
        endfor
    endfor
    savemat('kgrid.dat', kgrid, "%14.9f%14.9f%14.9f");
endfunction
