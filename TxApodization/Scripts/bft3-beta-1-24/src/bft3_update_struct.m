function opt = bft3_update_struct(opt,opt1)

  if nargin < 1, help(mfilename), error(mfilename), end
  if nargin == 1 && strcmp(opt, 'test')
    bft3_update_struct_test
    return
  end

  if strcmp(class(opt1),'struct')
    opt = bft3_update_struct_do(opt,opt1);
  end
end

% TODO: Support if b is a struct array
function opt = bft3_update_struct_do(opt,opt1)
  labels = fieldnames(opt1);
  
  for k=1:length(labels)
    eval(['b = opt1.',cell2mat(labels(k)),';']);
    if (strcmp(class(b),'struct'))
      eval(['check_exist = bft3_isvar(''opt.',cell2mat(labels(k)),''');']);
      if ~check_exist
        eval(['opt.',cell2mat(labels(k)),'=cell(0);']);
      end
      eval(['opt.',cell2mat(labels(k)),...
            ' = bft3_update_struct_do(opt.',...
            cell2mat(labels(k)),',opt1.',...
            cell2mat(labels(k)),');']);
    else
      eval(['opt.',cell2mat(labels(k)),'= opt1.', ...
            cell2mat(labels(k)),';']);
    end
  end
end

function bft3_update_struct_test
  opt.a = 2;
  opt.b = 3;
  opt.c.a = 1;
  opt.d.a.a = 1;
  
  opt1.a = 20;
  opt1.b = 3;
  opt1.c.a = 2;
  opt1.d.a.a = 2;
  
  opt = bft3_update_struct(opt,opt1);

end
