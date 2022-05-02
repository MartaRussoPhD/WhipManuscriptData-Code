function y = mydiff(t,x,varargin)
% Takes vector of time (may be not equidistantly spaced) Nx1 or 1xN,
% Takes vector of values or a matrix of values N x M
% Preferably, N is a column for both. Automatic reshaping works if M < N.



warning('off','MATLAB:interp1:NaNstrip');
s = size(x);

flag_filter = 0;
if any(strcmpi(varargin,'LPFilter'))
	flag_filter = 1;
    FS = varargin{find(strcmpi(varargin,'LPFilter'))+1};
    Fcut = varargin{find(strcmpi(varargin,'LPFilter'))+2};
end


if any(strcmpi(varargin,'Method'))
	InterpMethod = varargin{find(strcmpi(varargin,'Method'))+1};
else
	InterpMethod = 'spline';  %Could be changed (e.g. to makima) to avoid undulations
end

flag_Speed = 0;
if any(strcmpi(varargin,'Speed'))
	flag_Speed = 1;
end

if any(strcmpi(varargin,'DerivMethod'))
	DerivMethod = varargin{find(strcmpi(varargin,'DerivMethod'))+1};
else
	DerivMethod = 'numeric_spline';  % Could be changed (to 'analytic_spline') to ensure continuity. 
end


if any(strcmpi(varargin,'NoReshape'))
    ndof = s(2);
    %fprintf('<strong>Make sure mydiff arguments are time-columns!</strong>');
else
    
%Reshapes automatically.
ndof = min(s);
% columns expected
if s(1) < s(2) %I WANT IT TO BE COLUMN(S)
    x = x';
end
y = nan(size(x,1),size(x,2));
if flag_Speed
    y = nan(size(x,1),1);
end

if (nnz(isnan(x)) > 0.8 * numel (x))
    %disp(sprintf('<strong>More than 80%% nans</strong>'));
    return
end
st = size(t);
if st(1) < st(2)
    t = t';
end
assert(max(st) == max(s),'Time and data lengths are not equal!');
end



% Remember nan locations
ind0 = isnan(x(:,1));   
ind = isnan(diff(x,1,1));


dt = diff(t,1,1); 
yy = diff(x,1,1) ./ (dt * ones(1,ndof)); % Numerical derivatives - but these are shifted forward by one value
% Now we'll account for that shift by 

%t1 = t(1:end-1,1) + dt/2;  %works even if T not uniform

% To avoid edge effects with the next step
t1 = linspace(t(1),t(end-1),length(t)-3)';
yy = yy(2:end-1,:);

if strcmpi(DerivMethod,'numeric_spline')
    for idof = 1:ndof
        y(:,idof) = interp1My(t1,yy(:,idof),t,InterpMethod); % ,'extrap' % Extrap sometimes shoots the result in infinity.
    end
elseif strcmpi(DerivMethod,'analytic_spline') % To be continued maybe. Still not understand completely.
    for idof = 1:ndof
        ppspl = interp1My(t,x(:,idof),InterpMethod,'pp'); % ,'extrap' % I think this one is not correct
        %y(:,idof) = interp1(t1,yy(:,idof),t,InterpMethod); % ,'extrap' % I think this one is not correct
    end
end


if flag_filter
    y = LowFilt(y,Fcut - Fcut/10,Fcut + Fcut/10,'FS',FS,'FilterType','Butter4');
end

warning('on','MATLAB:interp1:NaNstrip');

y(ind0,:) = nan;

if flag_Speed
    y = sqrt(sum(y.^2,2));
end



% function v = velocity(data,time)
% ndof = size(data,1);
% dtime =  diff(time,1,2);
% data = diff(data,1,2)./(ones(ndof,1)*dtime);  % along rows
% timei = time(1:end-1) + dtime/2;
% 
% v.data = interp1(timei,data,time,'spline');
% v.time = time;
% end