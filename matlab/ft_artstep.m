% Altered by Berry van den Berg (berryv.dberg@gmail.com) to be used
% standalone without ERPlab


function trials = ft_artstep(data, twin, ampth, winms, stepms, chan)

fs       = data.fsample;
nch      = length(data.elec.label);
ntrial   = length(data.trial);
winpnts  = floor(winms*fs/1000);
stepnts  = floor(stepms*fs/1000);

tmp = abs(data.time{1}-twin(1)/1000);
[val p1] = min(tmp); %index of closest value

tmp = abs(data.time{1}-twin(2)/1000);
[val p2] = min(tmp); %index of closest value

trials = zeros(length(data.trial),1);

for ch=1:length(chan)
    fprintf('%g ',chan(ch));
    for i=1:ntrial;
        for j=p1:stepnts:p2-(winpnts-1)
            
            w1  = data.trial{i}(chan(ch), j:j+round(winpnts/2)-1);
            w2  = data.trial{i}(chan(ch), j+round(winpnts/2):j+winpnts-1);
            vs = abs(mean(w1)-mean(w2));
            
            if vs>ampth
       
                trials(i) = 1; % marks epoch with artifact

                break
            end
        end
    end
end

return