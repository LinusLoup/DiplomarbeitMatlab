function [x,y]=my_cirsg(bs,s)
%CIRSG  Gives geometry data for the cirsg PDE model
%
%   NE=CIRSG gives the number of boundary segment
%
%   D=CIRSG(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left hand region.
%   Row 4 contains the number of the right hand region.
%
%   [X,Y]=CIRSG(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.

% Copyright 1994-2003 The MathWorks, Inc.
% $Revision: 1.8.4.3 $

nbs=3;

if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 % start parameter value
  1 1 1 % end parameter value
  1 1 1 % left hand region
  0 0 0 % right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
  error(message('pde:cirsg:InvalidBs'))
end

if nargin==1,
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1,
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
  error(message('pde:cirsg:SizeBs'));
end

if ~isempty(s),

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=interp1([d(1,1),d(2,1)],[-0.70710678118654757 0],s(ii));
y(ii)=interp1([d(1,1),d(2,1)],[0.70710678118654757 0],s(ii));
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=interp1([d(1,2),d(2,2)],[0 -0.70710678118654757],s(ii));
y(ii)=interp1([d(1,2),d(2,2)],[0 -0.70710678118654757],s(ii));
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=1*cos(4.712388980384689674*s(ii)+(-2.356194490192344837))+(0);
y(ii)=1*sin(4.712388980384689674*s(ii)+(-2.356194490192344837))+(0);
end

end

