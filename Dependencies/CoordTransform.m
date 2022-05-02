function Coords2 = CoordTransform(Coords1,Rot,Trans,varargin)

% If a rotation matrix and translation vector are defined in the default coordinate frame,
% Perform translation and rotation into the new coordinate frame
Coords1 = squeeze(Coords1);


    
if size(Coords1,1) < size(Coords1,2)
    Coords1 = Coords1';
end

if any(strcmpi(varargin,'Reverse'))
    flag_reverse = 1;
else
    flag_reverse = 0;
end


Coords2 = Coords1;

if any(isnan(Rot),'all')
    s1 = max(size(Coords1));
    s2 = min(size(Coords1));
    Coords2 = nan(s1,s2);
    return
end

% If only 2 dimensions?
flag_2D = 0;
if size(Coords1,2) == 2
    flag_2D = 1;
    if numel(Rot) ~= 4 || numel(Trans) ~= 2
        disp('Coord error check arguments');
        return
    end
    Rot(1:3,3) = [0;0;1];
    Trans(end+1) = 0;
    Coords1(:,end+1) = zeros(length(Coords1),1);
    Coords2 = Coords1;
end


if ~flag_reverse
    % What if time column?
    if size(Coords1,2) == 3
        for iframe = 1:size(Coords1,1)
            Coords2(iframe,:) = Rot \ (Coords1(iframe,:)' - Trans);
        end
    else
        Coords2 = Rot \ (Coords1 - Trans);
    end
else
    if size(Coords1,2) == 3
        for iframe = 1:size(Coords1,1)
            Coords2(iframe,:) = Rot * Coords1(iframe,:)' + Trans;
        end
    else
        Coords2 = Rot * Coords1 + Trans;
    end
end

if flag_2D
    Coords2 = Coords2(:,1:2);
end

end