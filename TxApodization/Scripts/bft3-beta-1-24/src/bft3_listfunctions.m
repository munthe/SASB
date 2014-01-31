% function bft3_listfunctions(obj)
% 
% print list of non-inherited functions for an object
%
% $Id: bft3_listfunctions.m,v 1.3 2011-03-10 18:43:17 jmh Exp $

% @file  bft3_listfunctions.m
% @brief Function used for listing member functions
% ======================================================================
function bft3_listfunctions(obj)
  [m n] = methods(class(obj),'-full');
  for k=1:length(m)
    if (isempty(findstr(m{k},'L ')) && ...
        isempty(findstr(m{k},'TF ')) && ...
        isempty(findstr(m{k},'HM ')) && ...
        isempty(findstr(m{k},'Inherited')))
      fprintf('    ');
      fprintf(cell2mat(m(k)));
      fprintf('\n');
    end
  end
end
