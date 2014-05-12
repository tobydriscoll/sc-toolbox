function scgexprt(data)
%   Export data to base workspace.

%   Copyright 1997 by Toby Driscoll. Last updated 04/29/97.

tag = {'pol','map','phy','can'};
field = {'polygon','map','phypoints','canpoints'};

for j = 1:4
  name = get(findobj(gcf,'tag',tag{j}),'string');
  if ~isempty(name)
    assignin('base',name,getfield(data,field{j}))
  end
end
