function [data]  = contraMinusIpsi(layoutmw64, dataLeft, dataRight, powVSerp)
%-------- channel index for left and right ----%
chanIdx = zeros(length(layoutmw64.channels.conIpsi),2);
for iChans = 1:length(layoutmw64.channels.conIpsi)
    for i = 1:2
        chanIdx(iChans,i) = find(strcmp(dataLeft.label,layoutmw64.channels.conIpsi(iChans,i)));
    end
end
%-------- channel index for middle electrodes ----%
chanIdxMiddle =zeros(length(layoutmw64.channels.middle),1);
for i = 1:length(layoutmw64.channels.middle)
    chanIdxMiddle(i) = find(strcmp(dataLeft.label,layoutmw64.channels.middle(i)));
end;
%----- put it in a more convenient structure-----%
electrodes.left =  chanIdx(:,1);
electrodes.right =  chanIdx(:,2);
electrodes.middle = chanIdxMiddle;

switch powVSerp
    case 'power'
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
    case 'erp'
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