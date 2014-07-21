function data = scgimprt(data)
%   Import data from base workspace.

%   Copyright 1997 by Toby Driscoll. Last updated 05/23/97.

tag = {'pol','map','phy','can'};
field = {'polygon','map','phypoints','canpoints'};

for j = 1:4
  name = get(findobj(gcf,'tag',tag{j}),'string');
  if ~isempty(name)
    val = evalin('base',name);
    data = setfield(data,field{j},val);
    if j == 1
      % New polygon means map is not current
      data.iscurrent = 0;
      data.phypoints = [];
      data.canpoints = [];
    elseif j == 2
      % New map is current and defines the polygon, too
      data.iscurrent = 1;
      data.polygon = polygon(data.map);
      data.phypoints = [];
      data.canpoints = [];
    end
  end
end
