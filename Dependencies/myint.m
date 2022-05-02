function y = myint(t,x)

%disp(sprintf('<strong>Use cumtrapz instead?</strong>'));
% But how does it handle nans?...


    if size(t,1) < size(t,2)
        t = t';
    end
    if size(x,1) < size(x,2)
        x = x';
    end
    %Column vectors expected
    % X could be multi-column matrix?
    

%ind = isnan(diff(x,1,1));

s = size(x,1);
ndof = size(x,2);
dt = diff(t,1,1);  %replace with nans where x is nan?

for idof = 1:ndof
    ind0(idof,:) = isnan(x(:,idof));
end

SS = zeros(length(x),ndof);
for iframe = 2:s
    for idof = 1:ndof
        if isnan(x(iframe,idof)) && ~isnan(x(iframe-1,idof))
            dS = x(iframe-1,idof) .* dt(iframe,1);
            SS(iframe,idof) = SS(iframe-1,idof) + dS;
        elseif isnan(x(iframe,idof)) && isnan(x(iframe-1,idof))
            SS(iframe,idof) = SS(iframe-1,idof);
        elseif ~isnan(x(iframe,idof)) && isnan(x(iframe-1,idof))
            dS = x(iframe,idof) .* dt(iframe,1);
            SS(iframe,idof) = SS(iframe-1,idof) + dS;
        else
            dS = 1/2 * (x(iframe-1,idof)+x(iframe,idof)) .* dt(iframe-1,1);
            SS(iframe,idof) = SS(iframe-1,idof) + dS;
        end
    end
end
SS(1,:) = []; %dangerous moment!!!!!

%yy = x .* (ones(ndof,1)*dt);


t1 = t(1:end-1,1) + dt/2;  %s
y = interp1(t1,SS,t,'spline','extrap');

for idof = 1:ndof
    y(ind0(idof,:),idof) = nan;
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