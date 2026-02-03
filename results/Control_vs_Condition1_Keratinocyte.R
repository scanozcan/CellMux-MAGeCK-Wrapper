pdf(file='Control_vs_Condition1_Keratinocyte.pdf',width=4.5,height=4.5);
gstable=read.table('Control_vs_Condition1_Keratinocyte.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("TP53","MYC","JUNB","RELA","PIK3CB","MTOR","RB1","MDM2","CDKN1A","AKT2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Condition1_Rep1,Condition1_Rep2,Condition1_Rep3_vs_Control_Rep1,Control_Rep2,Control_Rep3 neg.'


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
targetmat=list(c(1600.704866399808,1610.0633966455498,1593.6831472751069,446.4767980121546,446.2870929101037,487.7017334572336),c(1267.2621345444857,1196.799167860779,1172.0315559584328,362.6906636808248,362.4721998513769,526.7178721338123),c(1057.8493137296878,1067.818612380182,1039.1413624865902,331.7012715308809,330.905551816272,261.8188253296728),c(1145.928525917886,1145.0314619195872,1147.0521962832745,292.67759252724784,235.11710260629852,234.09683205947215),c(1564.754167547482,1513.5473347212933,1640.6443434644048,674.8800957098891,673.7846597837906,457.92625920405516),c(1047.9628715452982,1047.6379812505647,1046.1355831956346,432.70373483440176,281.9228221066265,282.353635159451),c(1459.598373404429,1410.889341583675,1414.8309320009728,355.80413209194836,355.94116922342414,444.5786328146993),c(1072.2295932706181,1097.6508497022248,1113.080267125059,347.7698452382592,347.23312838615385,410.69619659556514),c(1557.564027777017,1555.6634344700597,1556.7136949558726,430.40822430477624,473.4997205265734,474.35410706787775),c(916.7428207343089,915.1477507909037,917.242087271817,232.99431875698556,250.35617407152156,250.52467992329474),c(1719.3421726124832,1717.1084835069978,1720.5782944249117,393.6800558307687,576.9077054691584,578.054896708258),c(960.782426828408,959.8961067739681,961.2057603000959,289.2343267328097,239.4711230249337,240.25727500840563))
targetgene="TP53"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(2514.751384720192,2511.172447520198,2513.922757707943,583.0596745248703,770.6616140984229,772.1088495996624),c(1700.468055715012,1698.6826898669124,1699.5956322977786,565.8433455526792,379.88828152591753,379.8939818508978),c(1538.6899108795458,1536.3602220852085,1537.729381602752,553.2180376397391,446.2870929101037,447.658854289166),c(1887.4116897471065,1885.5728825020635,1885.4420682809573,578.4686534656194,547.5180676433711,548.2794224550795),c(1482.0675601871326,1479.3280036754206,1482.7747903174036,521.0808902249825,414.7204448749988,413.77641807003187),c(1237.602807991317,1238.0378488647796,1234.979542339832,291.52983726243514,274.30328637401493,274.13971122753975))
targetgene="MYC"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1316.6943454664338,1315.2506984041847,1311.91597013932,453.363329601031,272.12627616469734,272.08623024456193),c(1528.8034686951562,1526.7086158927827,1527.7376377326887,614.0490666748142,470.23420521259703,471.273885593411),c(1421.8501396094869,1422.2957852656327,1418.827629548998,386.79352424189227,286.27684252526166,287.4873376168956),c(1337.3659973065212,1336.308748278568,1336.8953298144784,448.77230854178003,348.32163349081264,349.0917671062304),c(2061.7725791808866,2059.3017939657248,2060.297586007066,592.2417166433721,700.9972874002605,703.3172366699054),c(1202.5508766102992,1201.186261584609,1205.004310729642,356.9518873567611,304.781429304461,305.968666463696))
targetgene="JUNB"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1444.3193263921905,1442.47641639525,1441.8086404501437,527.9674218138589,414.7204448749988,463.0599616614997),c(1510.8281192689933,1509.1602409974635,1508.7533243795683,432.70373483440176,502.8893583523607,424.043822984921),c(1764.2805461778905,1762.734258234828,1758.5469211311524,539.4449744619862,710.7938333421895,622.2047378422812),c(1512.6256542116096,1610.9408153903157,1615.6649837892464,548.6270165804882,556.2261084806414,443.5518923232104),c(1320.2894153516663,1318.7603733832486,1318.9101908483644,327.11025047162997,313.4894701417314,433.28448740832124),c(1226.817598335619,1294.1926485298015,1296.9283543342249,386.79352424189227,353.76415901410655,369.62657693600863))
targetgene="RELA"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1403.874790183324,1402.1151541360155,1403.840013743903,449.9200638065928,403.8353938284109,404.5357536466317),c(1692.3791484732387,1690.7859211640186,1689.6038884277152,549.7747718453008,601.9433228763106,601.669928012503),c(1142.3334560326534,1142.3992056852892,1140.05797557423,331.7012715308809,293.89637825787315,294.674521057318),c(1369.7216262736144,1368.7732418349087,1368.8689101986813,404.0098532140833,386.41931215387024,387.08116529132013),c(1464.990978232278,1463.5344662696332,1464.7896513512894,511.89884810648056,424.5169908169279,425.0705634764099),c(1768.7743835344313,1767.121351958658,1762.5436186791778,500.4212954583532,500.71234814304313,502.0761003380784))
targetgene="PIK3CB"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(2029.4169502137936,2026.8373004093842,2027.3248312358571,698.9829562709566,556.2261084806414,558.5468273699686),c(1576.438144674488,1575.844065599677,1571.7013107609675,502.71680598797866,487.65028688713767,488.7284739487225),c(1318.49188040905,1317.0055358937168,1319.9093652353706,268.57473196618037,408.18941424704605,409.66945610407623),c(1866.740037907019,1865.3922513724463,1862.4610573798116,568.1388560823046,630.2444555974391,631.4454022656814),c(1701.3668231863203,1701.3149461012104,1695.5989347497532,456.80659539546923,526.8364706548541,527.7446126253012),c(1584.5270519162614,1581.9859968130388,1586.6889265660625,526.8196665490461,522.4824502362189,524.6643911508345))
targetgene="MTOR"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1390.3932781137019,1388.953872964526,1392.8490954868332,359.2473978863866,407.10090914238725,406.5892346296095),c(1112.6741294794847,1156.4379056015448,1155.045591379325,347.7698452382592,255.79869959481553,255.6583823807393),c(1550.3738880065519,1580.2311593235067,1580.6938802440245,595.6849824378104,339.6135926535423,339.8511026828302),c(1400.2797202980914,1387.1990354749942,1386.8540491647952,304.15514517537525,420.1629703982927,419.9368610189654),c(1560.2603301909414,1558.2956907043576,1560.7103925038978,600.2760034970613,454.995133747374,455.87277822107734),c(1377.8105335153878,1374.9151730482706,1379.8598284557509,414.33965059739796,445.1985878054449,444.5786328146993))
targetgene="RB1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1574.6406097318718,1572.3343906206132,1573.6996595349801,383.35025844745405,468.05719500327945,468.1936641189443),c(1047.9628715452982,1046.7605625057988,1048.1339319696472,376.4637268585776,252.53318428083915,253.6049013977615),c(1802.9275474441408,1801.3406830045305,1797.5147222243997,521.0808902249825,517.039924712925,516.4504672189232),c(1279.8448791427998,1279.27652986878,1279.9423897551173,359.2473978863866,301.5159139904847,301.8617044977404),c(1812.8139896285304,1811.8697079417223,1810.503989255482,412.0441400677725,710.7938333421895,710.5044201103277),c(1410.166162482481,1407.3796666046112,1409.835060065941,332.84902679569365,364.6492100606945,365.519614970053))
targetgene="MDM2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1579.1344470884123,1577.598903089209,1574.6988339219865,451.06781907140555,452.8181235380564,451.76581625512165),c(1116.2691993647172,1115.1992245975443,1113.080267125059,355.80413209194836,321.10900587434287,321.3697738360297),c(1114.471664422101,1113.4443871080123,1113.080267125059,421.2261821862744,334.1710671302484,334.7174002253856),c(1292.4276237411138,1291.5603922955036,1288.9349592381743,347.7698452382592,387.50781725852903,388.1079057828091),c(1218.7286910938458,1217.8572177351623,1217.9935777607243,356.9518873567611,315.66648035104896,315.20933088709626),c(1514.4231891542258,1512.6699159765274,1513.7491963146,391.3845453011432,500.71234814304313,502.0761003380784))
targetgene="CDKN1A"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1295.1239261550384,1295.0700672745675,1289.9341336251807,432.70373483440176,383.1537968398939,382.9742033253645),c(1029.9875221191353,1030.0896063552455,1026.1520954555078,206.5959476662926,241.64813323425128,242.31075599138345),c(1313.9980430525093,1312.6184421698867,1312.9151445263265,472.8751691028475,415.8089499796576,416.85663954449865),c(1260.0719947740206,1259.0958987391627,1259.9589020149906,294.97310305687336,412.5434346656812,413.77641807003187),c(1411.9636974250973,1409.1345040941433,1413.8317576139664,300.71187938093703,322.19751097900166,322.3965143275187),c(1553.0701904204761,1552.153759490996,1553.7161717948534,437.2947558936527,524.6594604455365,524.6643911508345))
targetgene="AKT2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetgenelist=c("DUSP6","SPRY4","CBLB","STK11","DUSP1","SPRY2","NF1","SOCS3","PTEN","CBL")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Condition1_Rep1,Condition1_Rep2,Condition1_Rep3_vs_Control_Rep1,Control_Rep2,Control_Rep3 pos.'


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
targetmat=list(c(1365.2277889170737,1337.1861670233338,1366.8705614246685,2061.368455603675,2061.6286682237474,1978.5289270991352),c(1756.1916389361172,1724.1278334651254,1756.5485723571398,2545.72117735465,2546.0134397969086,2716.7553404796636),c(1486.5613975436734,1544.2569907881023,1482.7747903174036,2438.9799377270656,2438.251434435688,2404.626231067034),c(1803.8263149154488,1841.701945263765,1803.5097685464377,2678.8607880729273,2678.8110625652807,2659.2578729562842),c(1088.4074077541648,1129.2379245137997,1084.1042099018753,1703.2688129821013,1704.5989938956643,1721.843804226907),c(1393.9883479989344,1330.1668170652063,1392.8490954868332,2094.6533582832444,2095.37232646817,2201.3316137522293))
targetgene="DUSP6"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1377.8105335153878,1397.7280604121856,1374.8639565207193,2165.814184701634,2166.125158270991,2133.566741313961),c(2087.836835848823,2046.1405127942353,2085.276945682225,3152.8837124405877,3151.2222779871954,3145.932865922029),c(976.0614738406465,1023.9476751418837,972.1966785571656,1509.2981732287487,1508.6680750570822,1674.613741618417),c(1518.9170265107666,1469.676397482995,1516.746719475619,2207.1333742348925,2208.576857352684,2235.2140499713632),c(1213.336086265997,1284.541042337376,1210.0001826646737,2007.4239581574766,2007.2034129908077,1829.6515558332428),c(1544.9812831787028,1574.966646854911,1541.7260791507774,2372.410132367927,2371.852623051502,2521.67464709677))
targetgene="SPRY4"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1154.9162006309675,1177.495955475928,1153.0472426053125,1827.226381581877,1827.6000707221076,1753.6727594630631),c(1440.724256506958,1384.5667792406962,1441.8086404501437,2101.539889872121,2100.8148519914635,2183.877025396918),c(1052.456708901839,1085.3669872755013,1048.1339319696472,1653.9153365951536,1652.3507488720425,1702.3357348886175),c(1205.2471790242237,1252.0765487810352,1202.006787568623,1961.513747564967,1961.4861985951386,1724.9240257013737),c(1683.3914737601574,1717.9859022517637,1682.6096677186708,2613.4387379786017,2612.4122511810947,2644.8835060754395),c(1367.02532385969,1350.3474481948235,1363.8730382636495,2165.814184701634,2165.036653166332,2161.288734584162))
targetgene="CBLB"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(2013.239135730247,2006.656669279767,2011.3380410437558,3064.506557050007,3061.9648594051746,3288.649794238988),c(1396.684650412859,1446.8635101190798,1391.849921099827,2367.819111308676,2367.498602632867,2079.1494952650487),c(1221.4249935077703,1206.450774053205,1217.9935777607243,1803.1235210208094,1802.5644533149555,2072.9890523161152),c(1439.8254890356498,1450.3731850981437,1438.8111172891247,2160.0754083775705,2160.682632747697,2245.4814548862523),c(1563.855400076174,1552.153759490996,1562.7087412779106,2440.127692991878,2439.339939540347,2329.6741751883433),c(1153.1186656883513,1203.818517818907,1150.0497194442935,1885.7619000873265,1884.2023361643646,1859.4270300864212))
targetgene="STK11"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1740.9125919238786,1606.553721666486,1742.560130939051,2583.5971010934704,2581.9341082506485,2341.99506108621),c(1708.5569629567854,1647.7924026704864,1706.589853006823,2660.4967038359237,2659.2179706814227,2514.4874636563477),c(1588.122121801494,1535.4828033404426,1586.6889265660625,2418.3203429604364,2417.5698374471713,2393.332085660656),c(1588.122121801494,1517.0570097003572,1590.685624114088,2287.476242771784,2286.9492248881165,2479.578286945725),c(1617.7814483546626,1690.7859211640186,1613.6666350152336,2499.8109667621407,2497.030710087263,2546.316418892504),c(1418.2550697242543,1486.3473536335484,1412.83258322696,2188.769289997889,2187.895260364167,2288.604555528787))
targetgene="DUSP1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1449.7119312200393,1445.108672629548,1448.802861159188,2101.539889872121,2104.08036730544,2139.7271842628943),c(1386.798208228469,1338.941004512866,1385.8548747777888,2270.259913799593,2269.533143213576,2104.8180075522714),c(1753.4953365221927,1619.7150028379754,1753.5510491961209,2624.916290626729,2625.474312437,2731.1297073605083),c(1239.4003429339332,1267.8700861868226,1236.9778911138449,1994.7986502445365,1990.875836420926,1875.854877950244),c(1363.4302539744574,1345.9603544709937,1361.874689489637,2040.708860837046,2039.8585661305715,1984.6893700480687),c(1260.0719947740206,1339.8184232576318,1253.9638556929524,2199.0990873812034,2195.5147960967784,2058.6146854352705))
targetgene="SPRY2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1363.4302539744574,1339.8184232576318,1363.8730382636495,1989.0598739204727,1990.875836420926,2004.197439386358),c(1154.0174331596593,1135.3798557271614,1153.0472426053125,1836.4084237003788,1834.1311013500604,1779.341271750286),c(1976.389669406613,1950.501869614745,1977.3661118855402,2841.8420356763363,2842.0868282640995,3050.44600021356),c(1318.49188040905,1309.108767190823,1316.9118420743516,2108.4264214609975,2107.3458826194164,2070.9355713331374),c(1034.481359475676,1055.5347499534585,1034.1454905515584,1556.356139086071,1553.2967843480926,1667.4265581779946),c(1149.5235958031187,1206.450774053205,1147.0521962832745,1860.5112842614462,1858.0782136525536,1856.3468086119544))
targetgene="NF1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1280.743646614108,1296.8249047640995,1282.9399129161363,1878.87536849845,1877.6713055364119,1889.2025043395997),c(2030.3157176851016,1932.0760759746597,2032.3207031708887,2960.060827952048,2961.822389776566,3088.43539839865),c(1038.9751968322166,1041.496050037203,1039.1413624865902,1542.583075908318,1542.4117333015047,1601.715166722704),c(1231.3114356921599,1232.7733363961838,1231.982019178813,1901.8304737947046,1902.706922943564,1835.8119987821763),c(1502.73921202722,1486.3473536335484,1502.7582780575303,2356.3415586605483,2353.348036272303,2237.267530954341),c(852.0315628001224,912.5154945566059,850.2974033423925,1330.2483519179618,1331.2417429976995,1303.9604241909194))
targetgene="SOCS3"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1187.2718295980608,1281.908786103078,1186.0199973765216,1908.7170053835812,1909.2379535715168,1682.8276655503282),c(1175.587852471055,1188.0249804131195,1173.0307303454392,1824.9308710522514,1825.42306051279,1827.598074850265),c(1205.2471790242237,1161.7024180701405,1204.0051363426358,1860.5112842614462,1860.2552238618712,1671.5335201439502),c(1452.4082336339638,1453.0054413324417,1451.800384320207,2068.2549871925517,2067.071193747041,2163.3422155671396),c(2095.9257430905964,2156.6952746347474,2089.27364323025,3366.3661916957567,3366.746288709636,3534.0407717048383),c(1283.4399490280323,1235.4055926304818,1286.9366104641615,2036.1178397777949,2034.4160406072774,1914.8710166268224))
targetgene="PTEN"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
targetmat=list(c(1333.7709274212887,1313.495860914653,1332.8986322664532,2030.3790634537313,2028.9735150839836,1997.0102559459356),c(1376.0129985727715,1410.889341583675,1370.867258972694,2117.608463579499,2116.053923456687,2039.106616096981),c(1438.0279540930335,1428.4377164789946,1437.8119429021185,2141.711324140567,2142.178045968498,2248.561676360719),c(1838.8782462964666,1780.2826331301474,1837.481697704653,2630.6550669507924,2632.0053430649527,2851.258344864711),c(1384.1019058145448,1366.1409856006107,1382.8573516167698,2160.0754083775705,2159.5941276430385,2077.096014282071),c(1519.8157939820749,1489.857028612612,1517.7458938626253,2127.9382609628137,2129.1159847125923,2345.075282560677))
targetgene="CBL"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition1_Rep1","Condition1_Rep2","Condition1_Rep3")

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
Sweave("Control_vs_Condition1_Keratinocyte_summary.Rnw");
library(tools);

texi2dvi("Control_vs_Condition1_Keratinocyte_summary.tex",pdf=TRUE);

