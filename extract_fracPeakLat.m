
%%%%-----correlational analysis CRS-----%%%
%%this is a very simple function for extracting the mean amplitude
function [lat,vals]  = extract_fracPeakLat(time_window,channels,data, direction,percentage)
chanIdx = zeros(1,length(channels));
for iElectrodes = 1: length(channels)
   chanIdx(iElectrodes) = find(strcmp(data.label,channels{iElectrodes}));
end

int_interest = data.time>=time_window(1) & data.time<=time_window(2);
temp = squeeze(mean(data.individual(:,chanIdx,int_interest),2));

tempTime = data.time(int_interest);
% find the peak 
switch direction
    case 'neg'
        [vals, idx] = min(temp,[],2);
    case 'pos'
        [vals, idx] = max(temp,[],2);

end



vals = vals*percentage;
lat = zeros(size(temp,1),1);
for i =1:size(temp,1)
    id= false;
    x = 0;
    while id ==false
        if idx(i)-x <=1
            disp('whoops widen window maybe?')
            lat(i) = tempTime(idx(i)-x);
            break
            
        end
        switch direction
            case 'neg'
                if  temp(i,idx(i)-x)>=vals(i)
                    id =true;
                    lat(i)= tempTime(idx(i)-x);
                end
            case 'pos'
                if  temp(i,idx(i)-x)<=vals(i)
                    id =true;
                    lat(i)= tempTime(idx(i)-x);
                end
        end
        x = x+1;
    end % end while
end % end loop

end

