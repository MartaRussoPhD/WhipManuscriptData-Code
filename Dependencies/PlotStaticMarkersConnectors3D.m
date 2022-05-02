function [curFig,sp3d] = PlotStaticMarkersConnectors3D(obj,objConn,varargin)

% If plot is provided
if any(strcmpi(varargin,'Subplot'))
    sp3d = varargin{find(strcmpi(varargin,'Subplot'))+1};
else
    sp3d = [];
end
% If position data was changed to target rf.
if any(strcmpi(varargin,'TargetRF'))
    TargetRF = 1;
else
    TargetRF = 0;
end

if any(strcmpi(varargin,'FullScreen'))
    flag_fullscreen = 1;
else
    flag_fullscreen = 0;
end

FigName = 'Whipper Sillhoette at ThrowOnset';
if any(strcmpi(varargin,'FigNameSpec'))
    if ~isempty(varargin{find(strcmpi(varargin,'FigNameSpec'))+1})
        FigName = varargin{find(strcmpi(varargin,'FigNameSpec'))+1};
    end
end


if isempty(sp3d)
    if flag_fullscreen
        curFig = Create_Reuse_Figure([],FigName,[1930 50 1600 900]);
    else
        curFig = Create_Reuse_Figure([],FigName,[2200 200 800 600]);
    end
    sp3d = subplot('Position',[0.1 0.1 0.85 0.85]); hold on;
    if flag_fullscreen
        sp3d.ZAxis.FontSize = 16;
        sp3d.YAxis.FontSize = 16;
        sp3d.XAxis.FontSize = 16;
        sp3d.XLabel.FontSize = 18;
        sp3d.YLabel.FontSize = 18;
        sp3d.ZLabel.FontSize = 18;
    end
else
    subplot(sp3d);
end
axis equal;

if TargetRF
    % If Y along the target
%     xlim([-2 2]); xlabel('X (m)','fontweight','bold');
%     ylim([-0.5 4.5]); ylabel('Y (m)','fontweight','bold');
%     zlim([-3 2]); zlabel('Z (m)','fontweight','bold');
%     view(36,30);
    % If X along the target
    xlim([-0.5 4.5]); xlabel('X (m)','fontweight','bold');
    ylim([-2 2]); ylabel('Y (m)','fontweight','bold');
    zlim([-2 2]); zlabel('Z (m)','fontweight','bold');
    view(-40,30)
else
    xlim([-2.5 2]); xlabel('X (m)','fontweight','bold');
    ylim([-2.5 2]); ylabel('Y (m)','fontweight','bold');
    zlim([0 4]); zlabel('Z (m)','fontweight','bold');
    view(30,30)
end



% Sphere sizes;
[sx,sy,sz] = sphere(15); %size 1
yratio = 1; zratio = 1;
scaleC = 0.01; % 1 cm radius.

% Create Markers
for imark = 1:size(obj,1)
    Marker(imark) = surface(sp3d,0,'EdgeColor','none','FaceColor',obj(imark,4:6));
    scaleCM = scaleC .* obj(imark,7); % Multiply by marker size from the list?...
    Marker(imark).XData = scaleCM*sx + obj(imark,1);
    Marker(imark).YData = scaleCM*yratio*sy + obj(imark,2);
    Marker(imark).ZData = scaleCM*zratio*sz + obj(imark,3);
end

% Create connectors
for icon = 1:size(objConn,1)
    Connect(icon,1) = line(sp3d,[0 0],[0 0],[0 0],'Color',objConn(icon,3:5),'LineWidth',objConn(icon,6));
    xyz1 = obj(round(objConn(icon,1)),1:3);
    xyz2 = obj(round(objConn(icon,2)),1:3);
    xyz = cat(1,xyz1,xyz2);
    Connect(icon,1).XData = xyz(:,1);
    Connect(icon,1).YData = xyz(:,2);
    Connect(icon,1).ZData = xyz(:,3);
end

end