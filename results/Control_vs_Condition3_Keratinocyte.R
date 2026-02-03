pdf(file='Control_vs_Condition3_Keratinocyte.pdf',width=4.5,height=4.5);
gstable=read.table('Control_vs_Condition3_Keratinocyte.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("FZD1","ADAM10","FZD7","PSEN1","HEY1","APC","NOTCH3","RBPJ","TGFBR1","TCF7")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Condition3_Rep1,Condition3_Rep2,Condition3_Rep3_vs_Control_Rep1,Control_Rep2,Control_Rep3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1688.806591487754,1634.5160328458453,1687.3817823207096,657.551342046378,500.67167621657813,392.26497032841417),c(1424.5654324151624,1355.5165167723194,1428.476639180269,305.83783350994327,345.22504150552624,484.4902768319598),c(1371.5374447101185,1360.7806585850274,1370.497881256927,415.06563119206584,279.4224783456331,341.8484694398092),c(1422.7678735099066,1411.6673627745383,1424.478104151073,408.5119633311385,318.52255210730874,340.6187986864286),c(1358.0557529207006,1210.7526169228483,1362.500811198535,253.4084906225244,537.8644293069525,329.5517619060031),c(1359.8533118259563,1258.1298932372206,1360.5015436839371,219.5478733410664,402.4446616445638,243.4748091693605))
targetgene="FZD1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1485.6824351938571,1409.0352918681845,1484.456129589013,440.188024658954,417.7032270149738,395.95398258855596),c(1306.8253241209125,1409.0352918681845,1301.5231520032962,356.08262044371963,278.4688180099825,381.19793354798867),c(1463.2129488781604,1457.2899251513413,1461.464553171136,476.23319789405446,430.1008113784319,527.5287532002811),c(1536.9128639936453,1544.148265061024,1536.4370849685608,494.80192350001533,503.53265722353,434.07377594335486),c(1128.8669925005956,1013.347298946297,1131.5854132624665,282.8999959966975,276.56149733868125,304.9583468383909),c(1382.3227981416528,1380.0825118982902,1384.492753859113,298.19188767219464,344.2713811698756,418.08805614940695))
targetgene="ADAM10"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1026.4061349010192,1067.7434310109468,1026.6238687460716,187.87181201325086,255.5809699543675,239.78579690921867),c(932.9330718277215,973.8662353509867,930.6590280453677,247.94710073841827,223.15651854224626,237.32645540245747),c(1287.9509556157273,1291.4694580510381,1285.5290118865123,318.94516923179793,473.01552648271,365.21221375404076),c(1411.0837406257444,1427.4597882126625,1408.483964034289,614.9525009503501,263.2102526395725,370.1308967675632),c(1433.5532269414412,1395.8749373364144,1434.4744417240631,347.3443966291498,365.25190855418936,359.06385998713773),c(1594.4347489618285,1545.9029789985934,1595.4154766492018,367.0054002119319,381.46413426025003,683.6969388796184))
targetgene="FZD7"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1478.4921995728341,1383.591939773429,1480.457594559817,555.969490202004,278.4688180099825,423.00673916292936),c(1698.6931654666605,1688.0348079417104,1700.3770211655967,480.6023098013394,363.3445878828881,485.7199475853404),c(1826.319847739817,1891.5816246997545,1824.3316070706724,389.94323772517765,535.9571086356513,504.1650088860495),c(1269.9753665631702,1287.9600301758996,1267.5356042551302,478.41775384769693,414.8422460080219,410.71003162912325),c(1325.6996926260974,1329.1958077087793,1323.5150946638742,543.9544324569705,283.2371196882356,419.3177269027875),c(1072.2438869850403,1077.394357667578,1069.6081203099286,212.99420548013904,407.2129633228169,311.10670060529395))
targetgene="PSEN1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1573.7628215513876,1542.3935511234547,1574.423167745923,586.5532735529982,581.7328047468812,541.0551314874677),c(1057.8634157429944,997.5548735081729,1057.6125152223406,301.4687216026583,230.78580122745126,252.08250444302476),c(1424.5654324151624,1416.0541476184617,1423.478470393774,453.29536038080875,390.04707728110566,473.4232400515343),c(1606.1188818459907,1696.808377629557,1601.4132791929958,370.28223414239557,638.9524248859187,529.9880947070424),c(1086.624358227086,1144.9508442639979,1087.6015279413105,296.0073317185522,464.4325834618544,302.4990053316297),c(1401.197166646838,1322.1769519585018,1403.4857952477942,570.1691039006799,362.3909275472375,463.58587402448944))
targetgene="HEY1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1530.6214078252503,1529.2331965916846,1528.4400149101689,467.4949740794847,460.6179421192519,430.384763683213),c(1251.9997775106128,1250.2336805181587,1250.5418303810472,369.18995616557436,327.10549512816436,320.94406663233883),c(1536.9128639936453,1534.4973384043926,1539.4359862404579,308.0223894635857,513.0692605800363,419.3177269027875),c(1088.4219171323418,1087.0452843242094,1089.6007954559084,301.4687216026583,323.29085378556186,375.04957978108564),c(1521.6336132989716,1520.459626903838,1518.443677337179,403.0505734470323,545.4937119921575,344.3078109465704),c(1349.067958394422,1345.865590115688,1349.5055723536482,550.5081003178979,337.59575882032124,456.2078495042058))
targetgene="APC"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1291.5460734262388,1240.5827538615272,1292.5264481876052,444.5571365662389,269.8858749891269,434.07377594335486),c(1199.8705692581968,1238.828039923958,1198.5608750014994,360.4517323510045,477.78382816096314,363.98254300066014),c(1142.3486842900136,1099.3282818871949,1142.5813845927555,395.40462760928375,317.5688917716581,284.0539440309206),c(1459.6178310676491,1486.2427051212355,1462.464186928435,460.9413062185573,318.52255210730874,480.80126457181797),c(1828.1174066450728,1735.4120842560826,1833.3283108863634,611.6756670198865,500.67167621657813,557.0408512814157),c(1234.9229679106834,1288.8373871446843,1233.5480565069643,430.357522867563,254.62730961871688,373.819909027705))
targetgene="NOTCH3"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1237.619306268567,1226.5450423609725,1235.5473240215622,457.66447228809363,376.6958325819969,404.5616778622202),c(1900.019762855302,1853.8552750420135,1902.3030401399944,446.7416925198814,431.0544717140825,640.6584625112971),c(2151.6780095911035,2080.213372988459,2151.211845707445,778.794197473534,704.7549880458118,438.9924589568773),c(1184.5913185635231,1218.6488296419104,1182.5667348847153,423.8038550066356,284.19078002388625,339.389127933048),c(1198.073010352941,1281.8185313944068,1195.5619737296024,419.4347430993507,380.5104739245994,295.12098081134604),c(1773.2918600347732,1769.6290060386848,1771.3510179338255,511.1860931523337,391.0007376167563,623.4430719639686))
targetgene="RBPJ"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1420.0715351520232,1479.2238493709583,1420.479569121877,326.5911150695465,374.78851191069566,531.2177654604229),c(1564.775027025109,1519.582269935053,1567.4257314448298,466.4026961026634,513.0692605800363,437.76278820349665),c(1172.907185679361,1156.3564848581987,1174.5696648263233,391.03551570199886,365.25190855418936,355.3748477269959),c(1382.3227981416528,1380.9598688670749,1381.493852587216,481.6945877781606,269.8858749891269,362.7528722472795),c(1212.453481594987,1240.5827538615272,1213.5553813609843,382.297291887429,505.43997789483126,402.102336355459),c(1636.677383235338,1559.940690499148,1640.3989957276567,415.06563119206584,541.679070649555,527.5287532002811))
targetgene="TGFBR1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(947.3135430697673,951.0549541625852,946.6531681621517,269.7926602748428,279.4224783456331,215.1923818416065),c(1473.099522857067,1560.8180474679327,1471.460890744126,430.357522867563,453.9423197696975,486.949618338721),c(1509.9494804148094,1454.6578542449874,1511.446241036086,376.83590200332293,426.2861700358294,448.82982498392215),c(1719.3650928771015,1620.4783213452906,1722.3689638261746,500.26331338412143,637.0451042146175,521.3803994333781),c(1432.6544474888133,1419.5635754936004,1432.475174209465,468.5872520563059,439.63741473493815,389.80562882165293),c(1648.3615161195003,1665.2235267533088,1647.3964320287498,482.7868657549818,547.4010326634587,498.0166551191465))
targetgene="TCF7"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("NOXA","BCLXL","MCL1","BAD","BCL2L1","HRK","BIK","BCL2L2","BIM","BAX")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Condition3_Rep1,Condition3_Rep2,Condition3_Rep3_vs_Control_Rep1,Control_Rep2,Control_Rep3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1346.3716200365384,1344.9882331469032,1314.5183908481831,2049.11348451662,2083.7478333966155,2073.224890199706),c(1396.7032693836986,1394.9975803676296,1512.445874793385,2128.8497768245693,2290.6921262328015,2277.350235260887),c(1679.8187969614755,1677.5065243162942,1604.4121804648928,2430.3184984272275,2536.7364928306624,2587.2272651128005),c(1260.9875720368914,1260.7619641435747,1305.521687032492,1961.7312463709216,2001.733044530662,1994.525961983347),c(1358.0557529207006,1355.5165167723194,1379.494585072618,2160.525838152385,2313.5799742884165,2148.2348061559233),c(1335.586266605004,1333.5825925527026,1271.5341392843263,1962.823524347743,2069.4429283618565,2251.5271494398944))
targetgene="NOXA"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1349.96673784705,1349.3750179908266,1388.491288888309,2244.631242367619,2157.179679241714,2258.905173960178),c(1182.7937596582674,1181.7998369529541,1149.5788208938484,1839.3961129669444,1780.4838466597168,1864.1808621250027),c(1440.743462562464,1439.742785775648,1468.4619894722289,2498.0397329901434,2319.30193630232,2244.1491249196106),c(1190.882774731918,1190.5734066408008,1167.5722285252305,1804.443217708665,1824.3522220996456,1838.3577763040098),c(1207.9595843318475,1207.2431890477096,1117.5905406602803,1895.1022897848268,1881.571842238683,1755.969835827509),c(1431.7556680361854,1430.9692160878012,1577.42206901782,2362.5972638643116,2339.328803350983,2233.082088139185))
targetgene="BCLXL"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(998.5439718695555,997.5548735081729,1005.6315598427926,1528.096889572895,1601.1957035573994,1583.8159303542238),c(1749.024814813821,1748.5724387878527,1865.3165911199314,2759.0941694504163,3011.659339984674,2683.141583876488),c(1596.2323078670843,1594.1576122817503,1528.4400149101689,2356.043596003384,2360.309330735297,2598.294301893226),c(1296.938750142006,1295.8562428949617,1308.520588304389,2006.514643420592,2002.6867048663125,2010.511681777295),c(1900.019762855302,1897.723123481247,1786.3455242933105,2739.433165867634,2796.1321041276324,2851.6064770896314),c(1628.5883681616872,1627.4971770955678,1714.3718937677827,2394.2733251921272,2484.2851743698784,2610.591009427032))
targetgene="MCL1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1384.1203570469086,1381.8372258358595,1332.5117984795652,2131.0343327782116,2070.396588697507,2111.344683554505),c(1513.544598225321,1512.5634141847759,1567.4257314448298,2434.6876103345126,2396.5484234900205,2320.3887116292085),c(1307.7241035735403,1306.3845265203777,1304.5220532751932,2316.72158883782,1961.6793104333358,1934.2720950676974),c(1333.7887076997483,1332.7052355839178,1385.492387616412,2184.5559536424516,2182.9285083042805,1996.9853034901084),c(1512.645818772693,1510.8087002472064,1510.4466072787868,2429.2262204504063,2306.904351938862,2340.0634436832984),c(882.6014224805612,882.6211105973807,882.6766076950158,1464.7447669172639,1312.23662185526,1383.3795975531848))
targetgene="BAD"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1520.7348338463437,1519.582269935053,1540.435619997757,2326.552090629211,2431.833855909094,2362.197517244149),c(1229.5302911949161,1228.2997562985418,1208.5572125744893,1990.1304737682735,1877.7572008960806,1856.802837604719),c(1043.4829445009486,1042.3000789161913,1062.6106840088355,1553.2192830397832,1603.1030242287006,1647.7588095300155),c(1442.5410214677197,1441.4974997132174,1445.470413054352,2297.060585255038,2142.8747742069545,2239.2304419060883),c(1859.574687487048,1857.3647029171523,1836.3272121582604,2809.338956384193,2829.5102158754044,2857.754830856534),c(1029.1024732589028,1028.2623674156364,1014.6282636584835,1622.0327955795203,1563.0492901313744,1569.0598813136567))
targetgene="BCL2L1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1373.3350036153743,1373.0636561480128,1465.4630882003319,2254.4617441590103,2093.284436753122,2137.1677693754978),c(1397.6020488363265,1396.752294305199,1456.466384384641,2153.9721702914576,2119.9869261513395,2177.746904237058),c(1498.2653475306472,1496.7709887466517,1452.467849355445,2182.3713976888093,2262.0823161632825,2171.598550470155),c(1259.1900131316356,1258.1298932372206,1322.515460906575,2019.6219791424467,2078.9795317183625,1967.4732054089739),c(1536.0140845410174,1535.3746953731772,1555.4301263572418,2288.3223614404683,2419.4362715456355,2365.886529504291),c(1385.9179159521643,1383.591939773429,1329.5128972076682,2236.9852965298705,2057.999004334049,2084.2919269801314))
targetgene="HRK"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2078.8768739282464,2075.826588144536,2055.2470050067413,3119.545901801421,3179.5035590591838,3066.7988589312376),c(1297.837529594634,1295.8562428949617,1263.5370692259341,1992.3150297219158,1943.5597640559738,1875.2478989054282),c(1589.0420722460613,1588.0161135002577,1632.4019256692648,2413.9343287749093,2486.1924950411794,2689.289937643391),c(725.3150182706852,725.5742131849244,759.7216555472389,1089.0011428907621,1154.882666472907,1131.29709311016),c(1654.6529722878954,1652.9405291903233,1599.4140116783979,2465.2713936855066,2700.76607056257,2549.107471758002),c(1450.6300365413704,1449.3937124322792,1411.482865306186,2122.296108963642,2200.094394345992,2217.0963683452374))
targetgene="BIK"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1700.4907243719163,1698.5630915671263,1699.3773874082976,2563.576411599417,2601.585395654905,2479.016238815307),c(1849.6881135081414,1848.5911332293053,1886.3089000232103,2849.753241526578,2913.4323254126593,2889.72627044443),c(1420.0715351520232,1416.9315045872465,1397.4879927040001,2173.633173874239,2202.955375352944,2175.2875627302965),c(1201.6681281634526,1200.2243332974324,1269.5348717697282,1924.593795159,1842.4717684770076,1914.5973630136075),c(1425.4642118677903,1422.1956463999545,1294.5257157022031,2085.15865775172,2070.396588697507,2062.157853419281),c(1668.1346640773133,1666.1008837220936,1629.4030243973677,2588.6988050663053,2552.9487185367234,2554.026154771524))
targetgene="BCL2L2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1012.0256636589735,1010.715228039943,954.6502382205437,1512.8049978973977,1520.1345750270964,1725.228066992994),c(1234.9229679106834,1233.56389811125,1241.5451265653562,1833.9347230828382,1887.2938042525868,1925.664399794033),c(1594.4347489618285,1594.1576122817503,1654.3938683298427,2476.194173453719,2473.7949106777214,2510.9876784032026),c(1594.4347489618285,1593.2802553129657,1592.4165753773048,2310.1679209768927,2465.211967656866,2356.049163477246),c(1843.3966573397465,1842.4496344478127,1846.3235497312505,2744.8945557517404,3100.349751200182,2771.677878119892),c(1228.6315117422885,1227.4223993297571,1224.5513526912732,1827.3810552219109,1770.9472433032106,1871.5588866452863))
targetgene="BIM"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1071.3451075324124,1070.3755019173007,1074.6062890964236,1725.7992033775367,1535.3931403975064,1678.5005783645308),c(1090.2194760375976,1089.6773552305635,1090.6004292132075,1702.8613658642912,1646.9713996686294,1658.825846310441),c(1710.3772983508227,1708.2140182237579,1697.3781198936997,2563.576411599417,2578.69754759929,2549.107471758002),c(1283.457058352588,1280.9411744256222,1259.5385341967383,1978.11541602324,1882.5255025743338,1842.0467885641517),c(1701.3895038245441,1699.4404485359112,1709.3737249812875,2563.576411599417,2745.5881063381494,2550.3371425113824),c(952.7062197855345,951.9323111313699,938.6560981037597,1551.0347270861407,1438.1197861611424,1441.1741229620734))
targetgene="BAX"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition3_Rep1","Condition3_Rep2","Condition3_Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("Control_vs_Condition3_Keratinocyte_summary.Rnw");
library(tools);

texi2dvi("Control_vs_Condition3_Keratinocyte_summary.tex",pdf=TRUE);

