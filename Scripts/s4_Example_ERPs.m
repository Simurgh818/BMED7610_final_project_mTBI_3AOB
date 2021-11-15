%% Example P3a, P3b

V = logical( double(IDENTITY.DEMO(:,2)<4) );
SIZEALL=sum( V  );

figure;  subplot(3,3,[1:6]); hold on
plot(tx2disp,squeeze(nanmean( MEGA_ERP(V,NovSite,:,1) ,1)),'k','Linewidth',2);
plot(tx2disp,squeeze(nanmean( MEGA_ERP(V,TargSite,:,2) ,1)),'g','Linewidth',2);
plot(tx2disp,squeeze(nanmean( MEGA_ERP(V,StdSite,:,3) ,1)),'c','Linewidth',2);
shadedErrorBar(tx2disp,squeeze(nanmean( MEGA_ERP(V,NovSite,:,1) ,1)),...
    squeeze(nanstd( MEGA_ERP(V,NovSite,:,1) ,1)) ./ sqrt(SIZEALL),'k');
shadedErrorBar(tx2disp,squeeze(nanmean( MEGA_ERP(V,TargSite,:,2) ,1)),...
    squeeze(nanstd( MEGA_ERP(V,TargSite,:,2) ,1)) ./ sqrt(SIZEALL),'g');
shadedErrorBar(tx2disp,squeeze(nanmean( MEGA_ERP(V,StdSite,:,3) ,1)),...
    squeeze(nanstd( MEGA_ERP(V,StdSite,:,3) ,1)) ./ sqrt(SIZEALL),'c');
legend({'Std','Targ','Nov'},'Location','NorthWest');
title(['N DataSets =',num2str(SIZEALL)])

plot([0 0],[-6 6],'k:');  plot([-500 1000],[0 0],'k:');
plot([NovT1 NovT1],[4 5],'m:'); plot([NovT2 NovT2],[4 5],'m:');
plot([TargT1 TargT1],[5.1 6],'r:'); plot([TargT2 TargT2],[5.1 6],'r:');

MAPLIMS=[-5 5];

% subplot(3,3,7);  topoplot( squeeze(nanmean(mean( MEGA_ERP(V,:,ERPWINS_tx2disp(1,1):ERPWINS_tx2disp(1,2),1)   ,3) ,1)) ,...
%     BV_Chanlocs_60,'maplimits',MAPLIMS,'emarker2',{StdSite,'d','k'}); title('Std')
% subplot(3,3,8);  topoplot( squeeze(nanmean(mean( MEGA_ERP(V,:,ERPWINS_tx2disp(2,1):ERPWINS_tx2disp(2,2),2) ,3) ,1)) ,...
%     BV_Chanlocs_60,'maplimits',MAPLIMS,'emarker2',{TargSite,'d','k'}); title('Targ') 
% subplot(3,3,9);  topoplot( squeeze(nanmean(mean( MEGA_ERP(V,:,ERPWINS_tx2disp(3,1):ERPWINS_tx2disp(3,2),3)   ,3) ,1)) ,...
%     BV_Chanlocs_60,'maplimits',MAPLIMS,'emarker2',{NovSite,'d','k'}); title('Nov') 
