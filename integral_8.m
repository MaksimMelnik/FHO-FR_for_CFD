function q = integral_8(f,varargin)
%Iterate INTEGRAL2 and INTEGRAL3 functions to calculate integrals of order
%4, 5, 6, 7, and 8. For convenience, the function is a straight pass-through to
%INTEGRAL, INTEGRAL2, and INTEGRAL3 for integrals of order 1, 2, and 3,
%respectively.
%
% Notes:
% *  This code is a straightforward extension of integralN.m (v.1.3.0.1)
%    of Mike Hosea (https://www.mathworks.com/matlabcentral/fileexchange/47919-integraln-m).
% *  This code has not been thoroughly tested.
% *  As with INTEGRAL, INTEGRAL2, and INTEGRAL3, your integrand MUST be
%    capable of accepting arrays and returning an array of the same size
%    with the value of the integrand computed "element-wise". In other
%    words, the method performs function evaluations in "batches". If your
%    integrand function only works with scalar arguments, instead of
%    passing f to this function directly, pass in...
%    ...for a 4D integral: @(x,y,z,w)arrayfun(f,x,y,z,w)
%    ...for a 5D integral: @(x,y,z,w,v)arrayfun(f,x,y,z,w,v)
%    ...for a 6D integral: @(x,y,z,w,v,u)arrayfun(f,x,y,z,w,v,u)
     ...for a 7D integral: @(x,y,z,w,v,u,s)arrayfun(f,x,y,z,w,v,u,s)
     ...for an 8D integral: @(x,y,z,w,v,u,s,t)arrayfun(f,x,y,z,w,v,u,s,t)
%    The same principle applies to functions that are used to define the
%    region of integration.
% *  Iterated adaptive quadrature is not a very efficient way of
%    approaching high order integrals. If speed is important, consider
%    sparse grid and Monte Carlo methods, among others, instead. Expect
%    this function to be slow on order 4 integrals, painfully slow on order
%    5, and excruciatingly slow on order 6. When using infinite limits (or
%    the 'iterated' method to deal with non-smoothness of the integrand) it
%    will take even longer. For example, the 4D hypersphere problem below
%    took less than 2 seconds on my machine, the 5D hypersphere problem
%    took over 1.5 minutes, and the 6D hypersphere problem took over 24
%    minutes. Do not expect 7D or 8D integrals to finish within your
%    lifetime on a commercial PC (year 2018).
%    The example integrating over all of R^4 took nearly 8
%    minutes.
% *  I recommend that you loosen the tolerances considerably, since for
%    technical reasons this method is likely to produce more accuracy than
%    requested on smooth problems. See the examples below. Always supply
%    both RelTol and AbsTol, never just one or the other. For an
%    explanation of these tolerances see
%
%    http://blogs.mathworks.com/loren/2013/12/26/double-integration-in-matlab-understanding-tolerances/
%
% Examples:
%    % Calculate the "volume" of the unit hypersphere.
%    xmin = -1;
%    xmax = 1;
%    ymin = @(x)-sqrt(1 - x.*x);
%    ymax = @(x) sqrt(1 - x.*x);
%    zmin = @(x,y)-sqrt(1 - x.*x - y.*y);
%    zmax = @(x,y) sqrt(1 - x.*x - y.*y);
%    wmin = @(x,y,z)-sqrt(1 - x.*x - y.*y - z.*z);
%    wmax = @(x,y,z) sqrt(1 - x.*x - y.*y - z.*z);
%
%    % 4D
%    f4 = @(x,y,z,w)ones(size(x));
%    q4 = integralN(f4,xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax,'AbsTol',1e-5,'RelTol',1e-3)
%    error = q4 - pi^2/2
%
%    % 5D
%    vmin = @(x,y,z,w)-sqrt(1 - x.*x - y.*y - z.*z  - w.*w);
%    vmax = @(x,y,z,w) sqrt(1 - x.*x - y.*y - z.*z  - w.*w);
%    f5 = @(x,y,z,w,v)ones(size(x));
%    q5 = integralN(f5,xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax,vmin,vmax,'AbsTol',1e-5,'RelTol',1e-3)
%    error = q5 - 2*pi^(5/2)/(gamma(5/2)*5)
%
%    % 6D
%    umin = @(x,y,z,w,v)-sqrt(1 - x.*x - y.*y - z.*z  - w.*w - v.*v);
%    umax = @(x,y,z,w,v) sqrt(1 - x.*x - y.*y - z.*z  - w.*w - v.*v);
%    f6 = @(x,y,z,w,v,u)ones(size(x));
%    q6 = integralN(f6,xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax,vmin,vmax,umin,umax,'AbsTol',1e-5,'RelTol',1e-3)
%    error = q6 - pi^3/6

%    % 7D
%    smin = @(x,y,z,w,v,u)-sqrt(1 - x.*x - y.*y - z.*z - w.*w - v.*v - u.*u);
%    smax = @(x,y,z,w,v,u) sqrt(1 - x.*x - y.*y - z.*z - w.*w - v.*v - u.*u);
%    f7 = @(x,y,z,w,v,u,s)ones(size(x));
%    q7 = integralN(f7,xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax,vmin,vmax,umin,umax,smin,smax,'AbsTol',1e-3,'RelTol',1e-3)
%    error = q7 - pi^(7/2)/gamma(7/2+1)
% 
%    % 8D
%    tmin = @(x,y,z,w,v,u,s)-sqrt(1 - x.*x - y.*y - z.*z - w.*w - v.*v - u.*u - s.*s);
%    tmax = @(x,y,z,w,v,u,s) sqrt(1 - x.*x - y.*y - z.*z - w.*w - v.*v - u.*u - s.*s);
%    f8 = @(x,y,z,w,v,u,s,t)ones(size(x));
%    q8 = integralN(f8,xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax,vmin,vmax,umin,umax,smin,smax,tmin,tmax,'AbsTol',1e-3,'RelTol',1e-3)
%    error = q8 - pi^(8/2)/gamma(8/2+1)
%
%    % Calculate a rapidly decaying 4D function over R^4.
%    f = @(x,y,z,w)exp(-x.*x - y.*y - z.*z - w.*w);
%    q = integralN(f,-inf,inf,-inf,inf,-inf,inf,-inf,inf,'AbsTol',1e-3,'RelTol',1e-3)
%
% Author:  Gabor Hannak
% Date:  2018-10-05

narginchk(3,inf);
if ~isa(f,'function_handle') || nargin(f) == 0
    error('integralN:BadIntegrand', ...
        ['First argument must be a function handle accepting one or more inputs.\n', ...
        'To integrate a constant c, integrate @(x,y)c*ones(size(x)) for 2D,\n', ...
        '@(x,y,z)c*ones(size(x)) for 3D, etc.']);
end
switch(nargin(f));
    case 1
        q = integral(f,varargin{:});
    case 2
        q = integral2(f,varargin{:});
    case 3
        q = integral3(f,varargin{:});
    case 4
        q = integral4(f,varargin{:});
    case 5
        q = integral5(f,varargin{:});
    case 6
        q = integral6(f,varargin{:});
    case 7
        q = integral7(f,varargin{:});
    case 8
        q = integral8(f,varargin{:});
    otherwise
        error('integralN:NTooLarge','N >= 8 is not supported.');
end

%--------------------------------------------------------------------------

function q = integral4(f,varargin)
% integral4(f) = integral2(integral2(f)).
narginchk(9,inf);
zmin = varargin{5};
zmax = varargin{6};
wmin = varargin{7};
wmax = varargin{8};
anyisinf = false;
if ~isa(zmin,'function_handle')
    anyisinf = anyisinf || isinf(zmin(1));
    zmin = @(x,y)zmin(1)*ones(size(x));
end
if ~isa(zmax,'function_handle')
    anyisinf = anyisinf || isinf(zmax(1));
    zmax = @(x,y)zmax(1)*ones(size(x));
end
if ~isa(wmin,'function_handle')
    anyisinf = anyisinf || isinf(wmin(1));
    wmin = @(x,y,z)wmin(1)*ones(size(x));
end
if ~isa(wmax,'function_handle')
    anyisinf = anyisinf || isinf(wmax(1));
    wmax = @(x,y,z)wmax(1)*ones(size(x));
end
if anyisinf
    method_override = {'method','iterated'};
else
    method_override = {};
end
inner = @(x,y)integral2( ...
    @(z,w)f(x*ones(size(z)),y*ones(size(z)),z,w), ...
    zmin(x,y), ...
    zmax(x,y), ...
    @(z)wmin(x*ones(size(z)),y*ones(size(z)),z), ...
    @(z)wmax(x*ones(size(z)),y*ones(size(z)),z), ...
    varargin{9:end},method_override{:});
q = integral2( ...
    @(xv,yv)arrayfun(inner,xv,yv), ...
    varargin{1:4},varargin{9:end});

%--------------------------------------------------------------------------

function q = integral5(f,varargin)
% integral5(4) = integral3(integral2(f)).
narginchk(11,inf);
wmin = varargin{7};
wmax = varargin{8};
vmin = varargin{9};
vmax = varargin{10};
anyisinf = false;
if ~isa(wmin,'function_handle')
    anyisinf = anyisinf || isinf(wmin(1));
    wmin = @(x,y,z)wmin(1)*ones(size(x));
end
if ~isa(wmax,'function_handle')
    anyisinf = anyisinf || isinf(wmax(1));
    wmax = @(x,y,z)wmax(1)*ones(size(x));
end
if ~isa(vmin,'function_handle')
    anyisinf = anyisinf || isinf(vmin(1));
    vmin = @(x,y,z,w)vmin(1)*ones(size(x));
end
if ~isa(vmax,'function_handle')
    anyisinf = anyisinf || isinf(vmax(1));
    vmax = @(x,y,z,w)vmax(1)*ones(size(x));
end
if anyisinf
    method_override = {'method','iterated'};
else
    method_override = {};
end
inner = @(x,y,z)integral2( ...
    @(w,v)f(x*ones(size(w)),y*ones(size(w)),z*ones(size(w)),w,v), ...
    wmin(x,y,z), ...
    wmax(x,y,z), ...
    @(w)vmin(x*ones(size(w)),y*ones(size(w)),z*ones(size(w)),w), ...
    @(w)vmax(x*ones(size(w)),y*ones(size(w)),z*ones(size(w)),w), ...
    varargin{11:end},method_override{:});
q = integral3( ...
    @(xv,yv,zv)arrayfun(inner,xv,yv,zv), ...
    varargin{1:6},varargin{11:end});

%--------------------------------------------------------------------------

function q = integral6(f,varargin)
% integral6(f) = integral4(integral2(f)).
narginchk(13,inf);
vmin = varargin{9};
vmax = varargin{10};
umin = varargin{11};
umax = varargin{12};
anyisinf = false;
if ~isa(vmin,'function_handle')
    anyisinf = anyisinf || isinf(vmin(1));
    vmin = @(x,y,z,w)vmin(1)*ones(size(x));
end
if ~isa(vmax,'function_handle')
    anyisinf = anyisinf || isinf(vmax(1));
    vmax = @(x,y,z,w)vmax(1)*ones(size(x));
end
if ~isa(umin,'function_handle')
    anyisinf = anyisinf || isinf(umin(1));
    umin = @(x,y,z,w,v)umin(1)*ones(size(x));
end
if ~isa(umax,'function_handle')
    anyisinf = anyisinf || isinf(umax(1));
    umax = @(x,y,z,w,v)umax(1)*ones(size(x));
end
if anyisinf
    method_override = {'method','iterated'};
else
    method_override = {};
end
inner = @(x,y,z,w)integral2( ...
    @(v,u)f(x*ones(size(v)),y*ones(size(v)),z*ones(size(v)),w*ones(size(v)),v,u), ...
    vmin(x,y,z,w), ...
    vmax(x,y,z,w), ...
    @(v)umin(x*ones(size(v)),y*ones(size(v)),z*ones(size(v)),w*ones(size(v)),v), ...
    @(v)umax(x*ones(size(v)),y*ones(size(v)),z*ones(size(v)),w*ones(size(v)),v), ...
    varargin{13:end},method_override{:});
q = integral4( ...
    @(xv,yv,zv,wv)arrayfun(inner,xv,yv,zv,wv), ...
    varargin{1:8},varargin{13:end});

%--------------------------------------------------------------------------

function q = integral7(f,varargin)
% integral6(f) = integral4(integral2(f)).
narginchk(15,inf);
vmin = varargin{9};
vmax = varargin{10};
umin = varargin{11};
umax = varargin{12};
smin = varargin{13};
smax = varargin{14};
anyisinf = false;
if ~isa(umin,'function_handle')
    anyisinf = anyisinf || isinf(umin(1));
    umin = @(x,y,z,w,v)umin(1)*ones(size(x));
end
if ~isa(umax,'function_handle')
    anyisinf = anyisinf || isinf(umax(1));
    umax = @(x,y,z,w,v)umax(1)*ones(size(x));
end
if ~isa(smin,'function_handle')
    anyisinf = anyisinf || isinf(smin(1));
    smin = @(x,y,z,w,v,u)smin(1)*ones(size(x));
end
if ~isa(smax,'function_handle')
    anyisinf = anyisinf || isinf(smax(1));
    smax = @(x,y,z,w,v,u)smax(1)*ones(size(x));
end
if anyisinf
    method_override = {'method','iterated'};
else
    method_override = {};
end
inner = @(x,y,z,w,v)integral2( ...
    @(u,s)f(x*ones(size(u)),y*ones(size(u)),z*ones(size(u)),w*ones(size(u)),v*ones(size(u)),u,s), ...
    umin(x,y,z,w,v), ...
    umax(x,y,z,w,v), ...
    @(u)smin(x*ones(size(u)),y*ones(size(u)),z*ones(size(u)),w*ones(size(u)),v*ones(size(u)),u), ...
    @(u)smax(x*ones(size(u)),y*ones(size(u)),z*ones(size(u)),w*ones(size(u)),v*ones(size(u)),u), ...
    varargin{15:end},method_override{:});
q = integral5( ...
    @(xv,yv,zv,wv,vv)arrayfun(inner,xv,yv,zv,wv,vv), ...
    varargin{1:10},varargin{15:end});

%--------------------------------------------------------------------------

function q = integral8(f,varargin)
% integral6(f) = integral4(integral2(f)).
narginchk(17,inf);
vmin = varargin{9};
vmax = varargin{10};
umin = varargin{11};
umax = varargin{12};
smin = varargin{13};
smax = varargin{14};
tmin = varargin{15};
tmax = varargin{16};
anyisinf = false;
if ~isa(smin,'function_handle')
    anyisinf = anyisinf || isinf(smin(1));
    smin = @(x,y,z,w,v,u)smin(1)*ones(size(x));
end
if ~isa(smax,'function_handle')
    anyisinf = anyisinf || isinf(smax(1));
    smax = @(x,y,z,w,v,u)smax(1)*ones(size(x));
end
if ~isa(tmin,'function_handle')
    anyisinf = anyisinf || isinf(tmin(1));
    tmin = @(x,y,z,w,v,u,s)tmin(1)*ones(size(x));
end
if ~isa(tmax,'function_handle')
    anyisinf = anyisinf || isinf(tmax(1));
    tmax = @(x,y,z,w,v,u,s)tmax(1)*ones(size(x));
end
if anyisinf
    method_override = {'method','iterated'};
else
    method_override = {};
end
inner = @(x,y,z,w,v,u)integral2( ...
    @(s,t)f(x*ones(size(s)),y*ones(size(s)),z*ones(size(s)),w*ones(size(s)),v*ones(size(s)),u*ones(size(s)),s,t), ...
    smin(x,y,z,w,v,u), ...
    smax(x,y,z,w,v,u), ...
    @(s)tmin(x*ones(size(s)),y*ones(size(s)),z*ones(size(s)),w*ones(size(s)),v*ones(size(s)),u*ones(size(s)),s), ...
    @(s)tmax(x*ones(size(s)),y*ones(size(s)),z*ones(size(s)),w*ones(size(s)),v*ones(size(s)),u*ones(size(s)),s), ...
    varargin{17:end},method_override{:});
q = integral6( ...
    @(xv,yv,zv,wv,vv,vu)arrayfun(inner,xv,yv,zv,wv,vv,vu), ...
    varargin{1:12},varargin{17:end});

%--------------------------------------------------------------------------
