function [data]  = contraMinusIpsi(electrodes, dataLeft, dataRight)
%--- get labels ---%
label = dataLeft.label;
dimord = dataLeft.dimord;

%----- put it in a more convenient structure-----%
[~,electrodes.left] =  ismember(electrodes.left,label);
[~,electrodes.right] =  ismember(electrodes.right,label);
[~,electrodes.middle] =  ismember(electrodes.middle,label);
electrodes.left = electrodes.left(electrodes.left ~= 0);
electrodes.right = electrodes.right(electrodes.right ~= 0);
electrodes.middle = electrodes.left(electrodes.middle ~= 0);

switch dimord
    case 'subj_chan_time_freq'
        conIpsi = dataLeft;
        conIpsi.powspctrm(:,electrodes.left,:,:) = (1/2 *(dataRight.powspctrm(:,electrodes.left,:,:) + dataLeft.powspctrm(:,electrodes.right,:,:))-...
            1/2 *(dataLeft.powspctrm(:,electrodes.left,:,:)+ dataRight.powspctrm(:,electrodes.right,:,:)));
        
        conIpsi.powspctrm(:,electrodes.middle,:,:) = dataLeft.powspctrm(:,electrodes.middle,:,:)- dataLeft.powspctrm(:,electrodes.middle,:,:);
        
        conIpsi.powspctrm(:,electrodes.right,:,:) =  (1/2 *(dataRight.powspctrm(:,electrodes.right,:,:) + dataLeft.powspctrm(:,electrodes.left,:,:))-...
            1/2 *(dataLeft.powspctrm(:,electrodes.right,:,:)+ dataRight.powspctrm(:,electrodes.left,:,:)));
        
        
        conIpsi.powspctrm(:,33,:,:) = zeros(size(conIpsi.powspctrm(:,33,:,:)));
        selchan = ft_channelselection({'all' '-RM'},dataLeft.label);
        conIpsi = ft_selectdata(conIpsi,'channel',selchan);
        conIpsi.powspctrm(:,strcmp(conIpsi.label,'LVEOG'),:) = zeros(size(conIpsi.powspctrm(:,strcmp(conIpsi.label,'LVEOG'),:)));
%        conIpsi.powspctrm(:,strcmp(conIpsi.label,'RVEOG'),:) = zeros(size(conIpsi.powspctrm(:,strcmp(conIpsi.label,'LVEOG'),:)))

        data = conIpsi;
    case 'subj_chan_time'
        conIpsi = dataLeft;
        conIpsi.individual(:,electrodes.right,:) = 1/2 *((dataRight.individual(:,electrodes.right,:) + dataLeft.individual(:,electrodes.left,:))-...
            (dataLeft.individual(:,electrodes.right,:)+ dataRight.individual(:,electrodes.left,:)));
        conIpsi.individual(:,electrodes.middle,:) = dataLeft.individual(:,electrodes.middle,:)- dataLeft.individual(:,electrodes.middle,:);
        conIpsi.individual(:,electrodes.left,:) = 1/2 *((dataRight.individual(:,electrodes.left,:)+  dataLeft.individual(:,electrodes.right,:))-...
            (dataLeft.individual(:,electrodes.left,:)+ dataRight.individual(:,electrodes.right,:)));
        selchan = ft_channelselection({'all' '-LM' '-RM' '-LVEOG' '-RVEOG'},dataLeft.label);

        %actichamp system -> LM, RM, LVEOG
        
        
        conIpsi = ft_selectdata(conIpsi,'channel',selchan);
        %conIpsi.individual(:,strcmp(conIpsi.label,'LVEOG'),:) = zeros(size(conIpsi.individual(:,strcmp(conIpsi.label,'LVEOG'),:)))
        
        
        data = conIpsi;
    otherwise
        disp('input power or ERP measure')
end

end