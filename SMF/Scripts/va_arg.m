function opt = va_arg(opt,varargs)
%function opt = va_arg(varargin)
% splits argument varargin and returns a struct with (name,value) pairs,
% e.g. opt = va_arg(opt,{'a',2,'b',3}) returns
% opt =
%     a: 2
%     b: 3
%
% $Id: va_arg.m,v 1.2 2009/08/05 13:13:01 jmh Exp $

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && strcmp(opt, 'test'), va_arg_test, clear, return, end

for k=1:1:length(varargs)/2,
    eval(['arg','=varargs{',int2str(2*k-1),'};']);
    eval(['opt.',arg,'=varargs{',int2str(2*k),'};']);
end;

function va_arg_test
args = {'a', 1, 'b', 2, 'c', 3};
opt.a = 0;
opt.b = 0;
opt.c = 0;
opt = va_arg(opt, args);

if ((opt.a ~= cell2mat(args(2))) || (opt.b ~= cell2mat(args(4))) || (opt.c ~= cell2mat(args(6))))
    error 'Should not get here'
end
