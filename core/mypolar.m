function hpol = mypolar(theta,rho,x,y,cf,vmag)
%POLAR  Polar coordinate plot.
%   POLAR(THETA, RHO) makes a plot using polar coordinates of
%   the angle THETA, in radians, versus the radius RHO.
%   POLAR(THETA,RHO,S) uses the linestyle specified in string S.
%   See PLOT for a description of legal linestyles.
%
%   See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.

%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 5.20 $  $Date: 2000/06/01 02:53:22 $

line_style = 'auto';

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
		'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
	  'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold on;
    maxrho = max(abs(rho(:)));
    hhh=plot([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho]);
    set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = maxrho+0.1; rticks = max(ticks-1,2);
    if rticks > 5   % see if we can reduce the number
        if rem(rticks,2) == 0
            rticks = rticks/2;
        elseif rem(rticks,3) == 0
            rticks = rticks/3;
        end
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~isstr(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% TOMOKO USE THIS SECTION TO ANNOTATE

contourf(x, y, cf, vmag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spokes
    th = (1:2)*pi/2;
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];

% draw radial circles
    c82 = cos(80*pi/180);
    s82 = sin(68*pi/180);
%    rinc = (rmax-rmin)/3;
    rinc = 10;
    i = 0;
    hhh = plot(xunit*i,yunit*i,'+','color',tc,...
                   'handlevisibility','off');
    set(hhh,'MarkerSize',0.0001)
    for i=0:rinc:rmax
        hhh = plot(xunit*i,yunit*i,'--','color',tc,'linewidth',0.0001,...
                   'handlevisibility','off');
        %text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
        %    ['  ' num2str(i)],'verticalalignment','bottom',...
        %    'handlevisibility','off')
    end
    i = 40;
    text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
            ['  50^o'],'verticalalignment','bottom' , ...
            'handlevisibility','off', 'FontSize',18,'FontName','times')

% annotate spokes in local time
    rt = 1.1*rmax;
    for i = 1:length(th);
        text(rt*cst(i),rt*snt(i),int2str(((24/360)*i*90) + 6),...
             'horizontalalignment','center',...
             'handlevisibility','off','FontSize',18,'FontName','times');
        if i == length(th);
%            loc = int2str(6);
             loc = '06';
        else
%            loc = int2str(00);
             loc = '00';
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off','FontSize',18,'FontName','times')
    end

    ii = 1;
    loc = 'LT';
%    text(-rt*cst(ii),(-rt*snt(ii)) - 7.0,loc,'horizontalalignment', ...
%            'center','handlevisibility','off','FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% END ANNOTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set view to 2-D
    view(2);
% set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

% transform data to Cartesian coordinates.
xx = rho.*cos(theta);
yy = rho.*sin(theta);

% plot data on top of grid
if strcmp(line_style,'auto')
    q = plot(xx,yy);
else
    q = plot(xx,yy,line_style);
end
if nargout > 0
    hpol = q;
end
if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end
set(get(gca,'xlabel'),'visible','on')
set(get(gca,'ylabel'),'visible','on')


