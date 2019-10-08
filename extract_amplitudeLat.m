
%%%%-----correlational analysis CRS-----%%%
%%this is a very simple function for extracting the mean amplitude
function [lat,vals]  = extract_amplitudeLat(time_window,channels,data, direction,amplitude)

chanIdx = zeros(1,length(channels));
for iElectrodes = 1: length(channels)
   chanIdx(iElectrodes) = find(strcmp(data.label,channels{iElectrodes}));
end

int_interest = data.time>=time_window(1) & data.time<=time_window(2);
temp = squeeze(mean(data.individual(:,chanIdx,int_interest),2));

tempTime = data.time(int_interest);
vals = amplitude;
lat = zeros(size(temp,1),1);

for i =1:size(temp,1)
    id= false;
    x = 1;
    while id ==false
        if  x >= length(tempTime)
            disp('whoops adjust window maybe?')
            lat(i) = tempTime(x);
            break
            
        end
        switch direction
            case 'neg'
                if  temp(i,x)<=vals
                    id =true;
                    lat(i)= tempTime(x);
                end
            case 'pos'
                if  temp(i,x)>=vals
                    id =true;
                    lat(i)= tempTime(x);
                end
        end
        x = x+1;
    end % end while
end % end loop

end

