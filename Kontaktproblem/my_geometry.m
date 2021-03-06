function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',4);
set(ax,'DataAspectRatio',[1 1.175 1]);
set(ax,'PlotBoxAspectRatio',[1.5 1 1]);
set(ax,'XLim',[-1.5 1.5]);
set(ax,'YLim',[-0.34999999999999998 2]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pdeellip(0,1,1,1,...
1.5707963267948966,'E1');
pderect([-1 1 0 -0.34999999999999998],'R1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','E1+R1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(9,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))
pdesetbd(8,...
'neu',...
2,...
char('0','0','0','0'),...
char('0','0'))
pdesetbd(7,...
'neu',...
2,...
char('0','0','0','0'),...
char('0','0'))
pdesetbd(6,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))
pdesetbd(5,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))
pdesetbd(4,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))
pdesetbd(3,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))
pdesetbd(2,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))
pdesetbd(1,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))

% PDE coefficients:
pdeseteq(1,...
char('2*((1E3)./(2*(1+(0.3))))+(2*((1E3)./(2*(1+(0.3)))).*(0.3)./(1-2*(0.3)))','0','(1E3)./(2*(1+(0.3)))','0','(1E3)./(2*(1+(0.3)))','2*((1E3)./(2*(1+(0.3)))).*(0.3)./(1-2*(0.3))','0','(1E3)./(2*(1+(0.3)))','0','2*((1E3)./(2*(1+(0.3))))+(2*((1E3)./(2*(1+(0.3)))).*(0.3)./(1-2*(0.3)))'),...
char('0.0','0.0','0.0','0.0'),...
char('0.0','0.0'),...
char('1.0','0','0','1.0'),...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1E3';...
'0.3';...
'0.0';...
'0.0';...
'1.0'])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');
