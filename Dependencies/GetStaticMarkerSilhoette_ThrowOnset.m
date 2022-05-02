function [Rot,Tr,obj,objConn,objLabels] = GetStaticMarkerSilhoette_ThrowOnset(sCur,iThrowOnset,iMinDist,iHSmax,flag_noWhipMarkers,flag_noWhipConn)

%






% Have the Target location for coordinate transform.
MarksTarg = sCur.Markers.Target(1:2,:,:);
MarksTarg = permute(MarksTarg,[3 2 1]) ./ 1000;

% Coordinate transform. Zero at target-1. X horizontal along the target rod.
% Y horizontal to the right from the target (looking from the person)
% Z vertical.
% See another commented alternative below!
TargMean = squeeze(mean(MarksTarg,1,'omitnan'));
Tr = TargMean(:,1); % Translation vector
TargX = TargMean(:,1) - TargMean(:,2);
TargX(3) = 0;
TargX = TargX ./ sqrt(sum(TargX.^2));
TargZ = [0; 0; 1];
TargY = cross(TargZ,TargX);
TargY = TargY ./ sqrt(sum(TargY.^2));
Rot = [TargX TargY TargZ]; % Rotation matrix

if nargin == 1 
    obj = nan;
    objConn = nan;
    objLabels = nan;
    return
end




% Have main body markers at ThrowOnset for the silhoette.
Head(1,:,:) = sCur.Markers.Head1;
Head(2,:,:) = sCur.Markers.Head2;
Head(3,:,:) = sCur.Markers.Head3;
Head(4,:,:) = sCur.Markers.Head4;
HeadFront = squeeze(mean(Head(1:2,:,:),1,'omitnan'));
HeadRear = squeeze(mean(Head(3:4,:,:),1,'omitnan'));
HeadR = squeeze(mean(Head([1 4],:,:),1,'omitnan'));
HeadL = squeeze(mean(Head(2:3,:,:),1,'omitnan'));
Neck = sCur.Markers.C7_Neck;
ShR = sCur.Markers.ShR;
ShL = sCur.Markers.ShL;
Chest = sCur.Markers.Sternum;
HipR = sCur.Markers.hipR;
HipL = sCur.Markers.hipL;
ElbR = sCur.Markers.Elb;
WrLat = sCur.Markers.WrL;
WrMed = sCur.Markers.WrM;
HandLat = sCur.Markers.HandL;
HandMed = sCur.Markers.HandM;
Hd1 = sCur.Markers.Hd1;
Hd0 = sCur.Markers.Hd0;
w10 = squeeze(sCur.Markers.Whip(10,:,:));
Hd = .5 * (Hd0 + w10);
Whip = sCur.Markers.Whip;

MarksBody = cat(3,HeadFront,HeadRear,HeadR,HeadL,Neck,Chest,ShR,ShL,HipR,HipL,ElbR,WrLat,WrMed,HandLat,HandMed,Hd1,Hd);
MarksBody = permute(MarksBody,[2 1 3]) ./ 1000;
MarksWhip = permute(Whip,[3 2 1]) ./ 1000;




% Add virtual markers: left elbow goes downward
% Left wrist - elbow a bit flexed
% Chest is mean of shoulders, C7 and Stern
% Pelvis is mean of Hips
% Feet are right under the Hips

% Create rotation matrix to position virtual left elbow and wrist properly.
ChestRotX = squeeze(mean(MarksBody(:,:,6) - MarksBody(:,:,5),1,'omitnan'));
ChestRotX(3) = 0;
ChestRotX = ChestRotX ./ sqrt(sum(ChestRotX.^2));
ChestRotY = squeeze(mean(MarksBody(:,:,8) - MarksBody(:,:,7),1,'omitnan'));
ChestRotY(3) = 0;
ChestRotY = ChestRotY ./ sqrt(sum(ChestRotY.^2));
ChestRotZ = cross(ChestRotX,ChestRotY);
ChestRot = [ChestRotX' ChestRotY' ChestRotZ']; % Rotation at ShR, X facing forward, Y facing laterally from the left shoulder
ChestTr = squeeze(mean(MarksBody(:,:,7),1,'omitnan')); % Translation at ShL

% Create markers in the local r.f.
ShL_t = CoordTransform(MarksBody(:,:,8),ChestRot,ChestTr');
UpArmL = mean(MarksBody(:,:,7) - MarksBody(:,:,11),1,'omitnan');
UpArmL = sqrt(sum(UpArmL.^2));
ElbL_t = ShL_t;
% Angle for shoulder flexion is 20deg
AShflex = 20 ./ 180 * pi;
ElbL_t(:,3) = ShL_t(:,3) - UpArmL * cos(AShflex);
ElbL_t(:,1) = ShL_t(:,1) + UpArmL * sin(AShflex);
LowArmL = mean(MarksBody(:,:,11) - MarksBody(:,:,12),1,'omitnan');
LowArmL = sqrt(sum(LowArmL.^2));
WrL_t = ElbL_t;
AElbFlex = 20 ./180 * pi;
WrL_t(:,3) = ElbL_t(:,3) - LowArmL * cos(AShflex + AElbFlex);
WrL_t(:,1) = ElbL_t(:,1) + LowArmL * sin(AShflex + AElbFlex);
% Convert back to the lab r.f.
ElbL = CoordTransform(ElbL_t,ChestRot,ChestTr','Reverse');
WrL = CoordTransform(WrL_t,ChestRot,ChestTr','Reverse');


Chest = squeeze(mean(MarksBody(:,:,5:8),3,'omitnan'));
Pelvis = squeeze(mean(MarksBody(:,:,9:10),3,'omitnan'));
FootR = squeeze(MarksBody(:,:,9)); FootR(:,3) = 0;
FootL = squeeze(MarksBody(:,:,10)); FootL(:,3) = 0;

% Target tripod top
Target5 = squeeze(MarksTarg(:,:,1)) + 3 * squeeze(MarksTarg(:,:,2) - MarksTarg(:,:,1));
TargetFloor = Target5; TargetFloor(:,3) = 0;

% Add virtual markers to the package.
MarksVirt = cat(3,ElbL,WrL,Chest,Pelvis,FootR,FootL,Target5,TargetFloor);


% Coordinate transform. Have zero at the target-1. Y horizontal along with the target rod
% X horizontal to the left from the target (looking from the person)
% Z vertical
% TargMean = squeeze(mean(MarksTarg,1,'omitnan'));
% Tr = TargMean(:,1); % Translation vector
% TargY = TargMean(:,1) - TargMean(:,2);
% TargY(3) = 0;
% TargY = TargY ./ sqrt(sum(TargY.^2));
% TargZ = [0; 0; 1];
% TargX = cross(TargY, TargZ);
% TargX = TargX ./ sqrt(sum(TargX.^2));
% Rot = [TargX TargY TargZ]; % Rotation matrix



% NO TRANSFORM TO TRY
%Tr = [0;0;0];
%Rot = eye(3);


% Transform into the target r.f.
TargMean(:,1) = CoordTransform(TargMean(:,1),Rot,Tr);
TargMean(:,2) = CoordTransform(TargMean(:,2),Rot,Tr);
for imark = 1:size(MarksBody,3)
    MarksBody(:,:,imark) = CoordTransform(MarksBody(:,:,imark),Rot,Tr);
end
for imark = 1:size(MarksVirt,3)
    MarksVirt(:,:,imark) = CoordTransform(MarksVirt(:,:,imark),Rot,Tr);
end
for imark = 1:size(MarksWhip,3)
    MarksWhip(:,:,imark) = CoordTransform(MarksWhip(:,:,imark),Rot,Tr);
end




% Prepare single frame, at ThrowOnset.
if length(iThrowOnset) == 1 % It's throw onset or HSmax indeed
    irange = iThrowOnset - 10 : iThrowOnset + 10;
elseif length(iThrowOnset) == 2 % For full-trial animation. [iframe sz]
    if iThrowOnset(1) > 2
        irange = (iThrowOnset(1)-2):iThrowOnset(1);
    elseif iThrowOnset(1) == 2
        irange = (iThrowOnset(1)-1):iThrowOnset(1);
    else
        irange = iThrowOnset(1);
    end
    if iThrowOnset(1) < iThrowOnset(2) - 1
        irange = [irange iThrowOnset(1) + (1:2)];
    elseif iThrowOnset(1) == iThrowOnset(2) - 1
        irange = [irange iThrowOnset(1) + 1];
    end
end
MarksBody = squeeze(mean(MarksBody(irange,:,:),1));
MarksWhip = squeeze(mean(MarksWhip(irange,:,:),1));
MarksVirt = squeeze(mean(MarksVirt(irange,:,:),1));

ColorsBody = repmat([0.6 0.1 0.1],4,1); % Head dark red
ColorsBody = cat(1,ColorsBody,repmat([0.1 0.4 0.1],6,1)); % Body dark green
ColorsBody = cat(1,ColorsBody,repmat([0.15 0.7 0.15],5,1)); % Arm brighter green
ColorsBody = cat(1,ColorsBody,repmat([0.4 0.2 0],2,1)); % Handle dark brown
ColorsWhip = [1 0 1; 0.85 0.325 0.098; 0.929 0.694 0.125; 0 1 0; 0 0.447 0.7410]; % w1-w5
ColorsWhip(6:10,:) = ColorModify(ColorsWhip, 1, .5, .8); % w6-w10
ColorsTarg = repmat([1 0 0],2,1); % Target bright red
ColorsVirt = repmat([0.5 0.8 0.5],size(MarksVirt,2),1); % Virtual markers
%ColorsVirt = cat(1, ColorsVirt,repmat([1 0 0],2,1)); % Virtual Targets


% Pack for plotting
MarkSize = 4; % roughly in cm;
objTarg = TargMean';
objTarg(:,4:6) = ColorsTarg;
objTarg(:,7) = ones(size(TargMean,2),1) .* MarkSize;

MarkSize = 3; % roughly in cm;
objBody(1:length(MarksBody),1:3) = MarksBody';
objBody(:,4:6) = ColorsBody;
objBody(:,7) = ones(length(MarksBody),1) .* MarkSize;

MarkSize = 3;
objWhip(1:length(MarksWhip),1:3) = MarksWhip';
objWhip(:,4:6) = ColorsWhip;
objWhip(:,7) = ones(length(MarksWhip),1) .* MarkSize;

MarkSize = 3;
objVirt(1:size(MarksVirt,2),1:3) = MarksVirt';
objVirt(:,4:6) = ColorsVirt;
objVirt(:,7) = ones(size(MarksVirt,2),1) .* MarkSize;

obj = cat(1,objTarg,objBody,objWhip,objVirt);
objLabels = {'Targ1','Targ2','HeadFront','HeadRear','HeadR','HeadL','Neck','Chest','ShR','ShL','HipR','HipL',...
    'ElbR','WrLat','WrMed','HandLat','HandMed','Hd1','Hd0',...
    'w1','w2','w3','w4','w5','w6','w7','w8','w9','w10','ElbL','WrL','Chest','Pelvis','FootR','FootL','Targ-5','TargFloor'}';

if flag_noWhipMarkers % Do not include whip markers to plot
    obj = cat(1,objTarg,objBody,objVirt);
    objLabels = {'Targ1','Targ2','HeadFront','HeadRear','HeadR','HeadL','Neck','Chest','ShR','ShL','HipR','HipL',...
        'ElbR','WrLat','WrMed','HandLat','HandMed','Hd1','Hd0',...
        'ElbL','WrL','Chest','Pelvis','FootR','FootL','Targ-5','TargFloor'}';
end

% Pack connectors (use marker indices)
objConn = [...
    1 2 1 0 0 1; %Target
    
    3 5 0 0 0 1;
    3 6 0 0 0 1;
    4 5 0 0 0 1;
    4 6 0 0 0 1; % Heads
    
    4 7 0 0 0 1; % Neck
    4 32 0 0 0 1; % To Chest
    
    %7 9 0 0 0 1;
    %7 10 0 0 0 1;
    %8 9 0 0 0 1;
    %8 10 0 0 0 1; % Chest area
    32 7 0 0 0 1;
    32 8 0 0 0 1;
    32 9 0 0 0 1;
    32 10 0 0 0 1;
    
    
    %9 11 0 0 0 1;
    %10 12 0 0 0 1; % Vertical torso and hips together
    32 33 0 0 0 1;
    33 11 0 0 0 1;
    33 12 0 0 0 1;
    
    11 34 0 0 0 1;
    12 35 0 0 0 1; %Legs virtual
    
    11 12 0 0 0 1;
    10 30 0 0 0 1;
    30 31 0 0 0 1; %Left arm virtual
    
    2 36 1 0 0 1; % Target virtual
    36 37 1 0 0 1;
    
    9 13 0 0 0 1;
    13 14 0 0 0 1;
    13 15 0 0 0 1;
    14 16 0 0 0 1;
    15 17 0 0 0 1;
    16 17 0 0 0 1; %Right arm
    18 19 .4 .2 0 2; %Handle
    
    20 21 0 0 0 1;
    21 22 0 0 0 1;
    22 23 0 0 0 1;
    23 24 0 0 0 1;
    24 25 0 0 0 1;
    25 26 0 0 0 1;
    26 27 0 0 0 1;
    27 28 0 0 0 1;
    28 29 0 0 0 1;
    29 19 0 0 0 1];  %Whip


if flag_noWhipConn
    
    objConn = [...
    1 2 1 0 0 1; %Target
    
    3 5 0 0 0 1;
    3 6 0 0 0 1;
    4 5 0 0 0 1;
    4 6 0 0 0 1; % Heads
    
    4 7 0 0 0 1; % Neck
    4 22 0 0 0 1; % To Chest
    
    %7 9 0 0 0 1;
    %7 10 0 0 0 1;
    %8 9 0 0 0 1;
    %8 10 0 0 0 1; % Chest area
    22 7 0 0 0 1;
    22 8 0 0 0 1;
    22 9 0 0 0 1;
    22 10 0 0 0 1;
    
    
    %9 11 0 0 0 1;
    %10 12 0 0 0 1; % Vertical torso and hips together
    22 23 0 0 0 1;
    23 11 0 0 0 1;
    23 12 0 0 0 1;
    
    11 24 0 0 0 1;
    12 25 0 0 0 1; %Legs virtual
    
    11 12 0 0 0 1;
    10 20 0 0 0 1;
    20 21 0 0 0 1; %Left arm virtual
    
    2 26 1 0 0 1; % Target virtual
    26 27 1 0 0 1;
    
    9 13 0 0 0 1;
    13 14 0 0 0 1;
    13 15 0 0 0 1;
    14 16 0 0 0 1;
    15 17 0 0 0 1;
    16 17 0 0 0 1; %Right arm
    18 19 .4 .2 0 2]; %Handle


end
end

