% CERR UTILITY FUNCTIONS (can be found at: https://github.com/adityaapte/CERR)
function [iV,jV,kV] = find3d(mask3M)

%   This is the 3D equivalent of the builtin command, find.
%   It returns the i,j,k indices of non-zero entries.

indV = find(mask3M(:));
[iV,jV,kV] = fastind2sub(size(mask3M),indV);
iV = iV';
jV = jV';
kV = kV';
end

function varargout = fastind2sub(siz,ndx)

%"fastind2sub"
%   FAST ind2sub is a faster version of ind2sub, based off the ind2sub that
%   ships with matlab.  Only one line is changed, and the original doc for
%   ind2sub is included below.
%
%   JRA 2/26/04

nout = max(nargout,1);
if length(siz)<=nout
  siz = [siz ones(1,nout-length(siz))];
else
  siz = [siz(1:nout-1) prod(siz(nout:end))];
end
n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1
  varargout{i} = floor(ndx/k(i)) + 1;
  ndx = ndx - (varargout{i}-1) * k(i);
end
end