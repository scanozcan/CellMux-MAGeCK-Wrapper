pdf(file='Control_vs_Condition2_Endothelial.pdf',width=4.5,height=4.5);
gstable=read.table('Control_vs_Condition2_Endothelial.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ANGPT2","NRP2","TIE1","ANGPT1","FLT4","NOTCH4","SOX17","PODXL","HGF","VEGFC")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Condition2_Rep1,Condition2_Rep2,Condition2_Rep3_vs_Control_Rep1,Control_Rep2,Control_Rep3 neg.'


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
targetmat=list(c(210.1949296070211,245.839664713231,225.19244939901287,70.99462650905413,74.0877072628529,71.95222116249103),c(210.1949296070211,230.81016556045705,257.52777546656347,72.19792526344487,92.15787976598776,76.64475732526218),c(195.79801662023883,179.28045417951782,162.83146341159392,109.50018664955806,69.57016413706918,61.78505947648686),c(227.47122519115982,206.118845523757,215.95378480828415,72.19792526344487,46.07893988299388,86.02982965080449),c(240.90834397882327,229.73662990668748,153.5927988208652,78.2144190353986,89.44735389051752,135.3014593599016),c(165.08460224843665,208.26591683129612,209.0247863652376,125.14307045663777,92.15787976598776,81.33729348803332),c(330.1692044968733,150.29499152773948,242.51494550662926,79.41771778978936,93.96489701630124,90.72236581357564),c(225.55163679292218,263.0162351735441,227.50211554669505,103.48369287760431,49.69297438362085,107.92833174373654),c(182.36089783257538,84.80931664779585,190.54745718378012,83.0276140529616,131.00875064772768,77.42684668572403),c(267.78258155415017,241.54552209815273,215.95378480828415,77.01112028100786,72.28069001253941,139.21190616221088),c(185.24028042993183,172.8392402569004,177.84429337152812,49.33524893002066,74.99121588800965,89.15818709265191),c(212.11451800525873,200.75116725490918,266.7664400572922,83.0276140529616,78.60525038863662,111.83877854604583))
targetgene="ANGPT2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(173.722750040506,212.5600594463744,170.91529492848156,140.78595426371749,57.82455201003153,71.17013180202916),c(220.7526657973281,192.16288202475263,174.37979415000484,58.96163896514665,97.57893151692821,89.94027645311378),c(225.55163679292218,243.69259340569187,222.8827832513307,101.07709536882282,77.70174176347987,109.49251046466026),c(268.74237575326896,178.20691852574825,174.37979415000484,101.07709536882282,87.64033664020404,89.94027645311378),c(217.87328319997164,237.25137948307446,189.39262410993905,64.97813273710038,113.84208676974959,47.7074509881734),c(224.59184259380336,222.2218803303005,255.21810931888126,57.7583402107559,58.72806063518827,90.72236581357564),c(224.59184259380336,148.14792022020035,217.10861788212523,77.01112028100786,83.12279351442032,89.15818709265191),c(182.36089783257538,234.03077252176575,228.65694862053616,73.40122401783562,91.25437114083101,127.48056575528301),c(187.15986882816946,192.16288202475263,255.21810931888126,87.8408090705246,70.47367276222593,89.15818709265191),c(297.5362017268335,187.86873940967436,251.75361009735798,64.97813273710038,98.48244014208495,91.5044551740375),c(212.11451800525873,191.08934637098307,241.36011243278816,78.2144190353986,126.49120752194398,109.49251046466026),c(227.47122519115982,210.41298813883526,218.26345095596633,120.32987543907478,53.30700888424782,76.64475732526218))
targetgene="NRP2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(143.00933566870384,217.92773771522224,237.89561321126487,96.26390035125982,121.07015577100351,118.0954934297407),c(341.6867348862991,181.42752548705695,227.50211554669505,104.68699163199506,67.7631468867557,78.2089360461859),c(133.41139367751566,265.1633064810832,203.25062099603213,87.8408090705246,122.877173021317,86.81191901126634),c(240.90834397882327,260.86916386600495,191.70229025762123,74.60452277222636,99.3859487672417,74.2984892438766),c(231.3104019876351,249.0602716745397,204.40545406987323,103.48369287760431,80.4122676389501,123.57011895297371),c(239.94854977970445,275.8986630187789,162.83146341159392,51.741846438802156,79.50875901379335,68.04177436018173),c(167.0041906466743,206.118845523757,185.92812488841577,113.11008291273029,78.60525038863662,75.08057860433846),c(235.14957878411036,188.94227506344393,227.50211554669505,61.36823647392814,140.94734552445186,49.27162970909711),c(161.24542545196138,265.1633064810832,249.4439439496758,68.58802900027263,69.57016413706918,75.86266796480032),c(245.70731497441736,273.75159171123977,190.54745718378012,107.09358914077656,99.3859487672417,126.69847639482116),c(170.84336744314956,217.92773771522224,177.84429337152812,128.75296671981002,73.18419863769616,76.64475732526218),c(168.92377904491192,237.25137948307446,175.53462722384594,140.78595426371749,75.89472451316638,96.19699133680865))
targetgene="TIE1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(219.79287159820927,219.0012733689918,169.76046185464048,79.41771778978936,47.885957133307365,96.19699133680865),c(225.55163679292218,202.8982385624483,162.83146341159392,101.07709536882282,65.95612963644221,89.15818709265191),c(201.55678181495173,231.88370121422662,271.38577235265655,93.85730284247833,131.91225927288443,114.96713598789327),c(261.0640221603184,209.3394524850657,160.52179726391174,85.4342115617431,123.78068164647375,51.617897790482694),c(150.68768926165438,186.7952037559048,205.56028714371433,104.68699163199506,65.05262101128547,62.56714883694872),c(204.4361644123082,182.50106114082652,247.13427780199362,62.571535228318886,74.99121588800965,130.60892319713045),c(164.12480804931783,206.118845523757,192.8571233314623,145.5991492812805,97.57893151692821,116.53131470881698),c(156.4464544563673,170.69216894936127,165.14112955927612,128.75296671981002,102.09647464271193,65.69550627879615),c(242.8279323770609,207.19238117752656,204.40545406987323,74.60452277222636,89.44735389051752,91.5044551740375),c(176.60213263786247,264.08977082731366,249.4439439496758,74.60452277222636,128.29822477225747,77.42684668572403),c(140.1299530713474,228.66309425291792,189.39262410993905,78.2144190353986,78.60525038863662,71.17013180202916),c(241.8681381779421,147.07438456643078,310.6500968632537,115.5166804215118,119.26313852069003,56.31043395325384))
targetgene="ANGPT1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(161.24542545196138,240.47198644438316,218.26345095596633,74.60452277222636,147.27190590054906,78.99102540664775),c(154.52686605812966,239.3984507906136,230.96661476821834,85.4342115617431,85.83331938989056,135.3014593599016),c(228.43101939027864,202.8982385624483,200.94095484834995,79.41771778978936,93.06138839114449,92.28654453449936),c(177.56192683698129,199.6776316011396,257.52777546656347,108.2968878951673,83.12279351442032,98.54325941819422),c(139.17015887222857,215.7806664076831,267.92127313113326,96.26390035125982,130.10524202257093,94.63281261588493),c(205.395958611427,280.1928056338572,228.65694862053616,137.17605800054525,62.34209513581524,72.73431052295288),c(186.20007462905065,198.60409594737004,233.27628091590051,75.80782152661712,74.99121588800965,189.26562523176986),c(252.42587436824908,198.60409594737004,161.67663033775284,113.11008291273029,140.04383689929512,102.45370622050352),c(143.96912986782266,187.86873940967436,150.12829959934191,57.7583402107559,112.03506951943609,101.67161686004167),c(237.069167182348,216.85420206145267,206.7151202175554,113.11008291273029,84.92981076473382,112.62086790650768),c(217.87328319997164,198.60409594737004,289.86310153411404,68.58802900027263,102.99998326786867,134.51936999943973),c(159.32583705372375,288.7810908640137,262.1471077619278,78.2144190353986,79.50875901379335,53.96416587186827))
targetgene="FLT4"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(234.18978458499154,390.7669779721227,221.72795017748962,89.04410782491534,93.96489701630124,82.11938284849519),c(268.74237575326896,266.2368421348528,232.12144784205944,66.18143149149113,91.25437114083101,51.617897790482694),c(208.27534120878346,255.5014855971571,185.92812488841577,72.19792526344487,36.140345006269705,77.42684668572403),c(222.67225419556573,264.08977082731366,226.34728247285398,43.31875515806692,82.21928488926359,116.53131470881698),c(373.35994345722014,282.3398769413963,176.68946029768702,62.571535228318886,58.72806063518827,78.2089360461859),c(219.79287159820927,301.66351870924854,212.48928558676087,70.99462650905413,48.789465758464104,93.06863389496121))
targetgene="NOTCH4"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(224.59184259380336,187.86873940967436,200.94095484834995,55.3517427019744,29.81578463017251,22.68059145339391),c(251.46608016913027,146.0008489126612,329.1274260447111,69.79132775466337,53.30700888424782,58.65670203463942),c(318.6516741074475,269.4574490961615,158.21213111622956,60.16493771953739,40.65788813205342,75.86266796480032),c(211.1547238061399,199.6776316011396,182.46362566689248,46.92865142123917,61.4385865106585,67.25968499971987),c(174.68254423962483,170.69216894936127,271.38577235265655,83.0276140529616,109.32454364396587,55.52834459279199),c(211.1547238061399,283.4134125951659,236.7407801374238,57.7583402107559,56.017534759718046,50.83580843002083))
targetgene="SOX17"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(188.11966302728828,243.69259340569187,256.37294239272234,97.46719910565058,58.72806063518827,58.65670203463942),c(266.8227873550313,214.70713075391353,220.5731171036485,72.19792526344487,115.64910402006306,87.59400837172821),c(142.04954146958502,234.03077252176575,222.8827832513307,60.16493771953739,91.25437114083101,67.25968499971987),c(217.87328319997164,208.26591683129612,173.22496107616377,81.82431529857085,65.95612963644221,81.33729348803332),c(236.10937298322918,366.0756579354226,262.1471077619278,74.60452277222636,48.789465758464104,112.62086790650768),c(298.49599592595234,285.560483902705,167.4507957069583,69.79132775466337,84.02630213957707,104.7999743018891))
targetgene="PODXL"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(195.79801662023883,198.60409594737004,205.56028714371433,68.58802900027263,56.92104338487479,49.27162970909711),c(127.65262848280277,154.58913414281776,352.22408752153297,85.4342115617431,97.57893151692821,35.976110581245514),c(228.43101939027864,171.76570460313084,220.5731171036485,43.31875515806692,65.95612963644221,61.002970116025),c(248.5866975717738,263.0162351735441,169.76046185464048,39.70885889489468,44.27192263268039,84.46565092988077),c(265.86299315591253,275.8986630187789,256.37294239272234,75.80782152661712,52.403500259091075,74.2984892438766),c(246.66710917353618,127.75074279857856,244.82461165431144,22.862676333424208,64.14911238612873,46.92536162771154))
targetgene="HGF"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(180.44130943433774,215.7806664076831,222.8827832513307,33.69236512294094,57.82455201003153,81.33729348803332),c(206.35575281054582,185.72166810213523,204.40545406987323,80.62101654418011,45.17543125783713,81.33729348803332),c(189.0794572264071,207.19238117752656,190.54745718378012,68.58802900027263,71.37718138738268,53.96416587186827),c(265.86299315591253,228.66309425291792,194.0119564053034,104.68699163199506,96.67542289177146,82.90147220895705),c(180.44130943433774,222.2218803303005,215.95378480828415,114.31338166712105,66.85963826159896,53.96416587186827),c(251.46608016913027,336.0166596298747,162.83146341159392,66.18143149149113,81.31577626410684,107.92833174373654))
targetgene="VEGFC"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetgenelist=c("CCNE1","MDM2","MDM4","CDKN1A","TP53","CCND2","CCND1","CDKN2A","CDK2","CDK4")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Condition2_Rep1,Condition2_Rep2,Condition2_Rep3_vs_Control_Rep1,Control_Rep2,Control_Rep3 pos.'


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
targetmat=list(c(218.83307739909046,164.25095502674387,250.5987770235169,453.64363040531197,347.8508206853459,508.3580843002083),c(203.47637021318937,196.4570246398309,235.5859470635827,282.7752072818258,264.7280271709256,240.1014336617907),c(201.55678181495173,153.5155984890482,162.83146341159392,277.96201226426274,526.745528466381,282.33425912673107),c(143.00933566870384,268.38391344239193,181.3087925930514,406.7149789840728,353.2718724362864,391.8267695913913),c(230.35060778851627,333.8695883223356,264.45677390961,411.52817400163576,398.44730369412355,441.88048866095033),c(183.3206920316942,237.25137948307446,209.0247863652376,370.61601635235036,485.1841317091708,436.4058631377173))
targetgene="CCNE1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(167.9639848457931,207.19238117752656,213.64411866060198,371.8193151067411,322.55257918095714,251.0506847082567),c(174.68254423962483,196.4570246398309,213.64411866060198,427.1710578087155,379.47362256583193,324.56708459167146),c(189.0794572264071,208.26591683129612,247.13427780199362,380.2424063874763,270.1490789218661,467.68943755619165),c(207.31554700966464,259.7956282122354,251.75361009735798,448.83043538774893,318.93854468033015,386.35214406815834),c(161.24542545196138,181.42752548705695,244.82461165431144,335.7203524750187,392.1227433180263,401.2118419169336),c(168.92377904491192,263.0162351735441,206.7151202175554,346.55004126453537,308.999949803606,311.27156546381985))
targetgene="MDM2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(188.11966302728828,171.76570460313084,203.25062099603213,370.61601635235036,285.5087255495307,375.4028930216923),c(231.3104019876351,252.2808786358484,192.8571233314623,466.87991670361015,394.83326919349656,341.7730505218324),c(188.11966302728828,205.04530986998745,147.81863345165974,316.46757240476666,340.622751684092,364.4536419752263),c(290.8176423330018,188.94227506344393,155.90246496854738,315.2642736503759,355.98239831175664,395.73721639370063),c(190.03925142552592,269.4574490961615,267.92127313113326,291.198298562561,449.0437867029011,315.18201226612916),c(240.90834397882327,241.54552209815273,232.12144784205944,317.67087115915746,330.6841568073678,417.6357184866327))
targetgene="MDM4"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(168.92377904491192,244.76612905946143,175.53462722384594,237.04985461497733,313.51749292938973,301.1044037778157),c(171.80316164226838,215.7806664076831,230.96661476821834,203.3574894920364,315.3245101797032,292.50142081273526),c(324.4104393021604,179.28045417951782,128.18647119636117,255.09933593083855,386.70169156708585,301.88649313827756),c(250.50628597001145,193.2364176785222,244.82461165431144,373.02261386113184,381.2806398161454,307.36111866151055),c(205.395958611427,123.45660018350029,289.86310153411404,435.5941490894507,341.52626030924876,317.52828034751474),c(200.59698761583292,193.2364176785222,205.56028714371433,202.15419073764565,323.4560878061139,416.07153976570896))
targetgene="CDKN1A"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(285.0588771382889,286.6340195564746,194.0119564053034,395.88529019455603,203.2894406602671,319.0924590684385),c(232.2701961867539,179.28045417951782,190.54745718378012,339.3302487381909,341.52626030924876,295.6297782545827),c(210.1949296070211,201.82470290867874,182.46362566689248,341.7368462469724,356.8859069369134,310.489476103358),c(131.49180527927803,172.8392402569004,135.11546963940773,229.83006208863284,355.0788896865999,201.7790549991596),c(238.02896138146681,126.67720714480899,236.7407801374238,492.1491905458159,362.3069586878538,462.9969013934205),c(207.31554700966464,246.91320036700057,270.23093927881547,302.02798735207773,436.3946659507067,374.6208036612304))
targetgene="TP53"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(223.63204839468455,221.14834467653094,222.8827832513307,440.4073441070137,329.78064818221105,462.9969013934205),c(172.7629558413872,176.05984721820911,138.579968860931,309.2477798784222,389.4122174425561,345.6834973241417),c(190.99904562464474,309.1782682856355,227.50211554669505,391.07209517699306,385.7981829419291,320.65663778936215),c(229.39081358939745,242.6190577519223,213.64411866060198,388.6654976682116,293.64030317594137,320.65663778936215),c(197.71760501847646,200.75116725490918,299.10176612484275,269.53892098352753,252.98241504388795,319.0924590684385),c(226.511430992041,191.08934637098307,182.46362566689248,335.7203524750187,312.613984304233,455.9580971492638))
targetgene="CCND2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(186.20007462905065,224.36895163783964,356.8434198168973,412.7314727560265,291.8332859256279,423.8924333703275),c(277.38052354533835,163.1774193729743,244.82461165431144,416.34136901919874,435.49115732555,409.03273552155224),c(275.4609351471007,216.85420206145267,333.7467583400755,203.3574894920364,376.7630966903617,563.8864288930004),c(231.3104019876351,172.8392402569004,196.3216225529856,321.2807674223297,538.4911405934187,291.71933145227337),c(187.15986882816946,193.2364176785222,227.50211554669505,293.60489607134247,289.12276005015764,273.7312761616506),c(154.52686605812966,162.10388371920473,209.0247863652376,267.13232347474604,494.21921796073826,307.36111866151055))
targetgene="CCND1"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(236.10937298322918,258.7220925584658,180.1539595192103,336.9236512294094,481.57009720854387,401.2118419169336),c(193.8784282220012,244.76612905946143,241.36011243278816,435.5941490894507,294.5438118010981,450.48347162603073),c(300.41558432419,179.28045417951782,200.94095484834995,362.1929250716151,253.8859236690447,326.91335267305703),c(194.83822242112,224.36895163783964,180.1539595192103,406.7149789840728,362.3069586878538,333.95215691721376),c(229.39081358939745,215.7806664076831,196.3216225529856,273.14881724669976,279.1841651734335,328.47753139398077),c(193.8784282220012,198.60409594737004,160.52179726391174,452.44033165092117,303.5788980526655,279.20590168488366))
targetgene="CDKN2A"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(188.11966302728828,271.60452040370063,173.22496107616377,512.6052693704586,242.14031154200703,304.23276121966313),c(284.09908293917005,242.6190577519223,168.60562878079938,299.62138984329624,413.8069503217881,344.11931860321795),c(199.6371934167141,320.98716047710076,181.3087925930514,367.00612008917807,429.16659694945275,372.27453557984484),c(246.66710917353618,221.14834467653094,136.2703027132488,311.6543773872037,224.07013903887218,365.2357313356881),c(381.0382970501707,238.32491513684403,213.64411866060198,316.46757240476666,328.87713955705436,463.77899075388234),c(238.02896138146681,313.4724109007138,221.72795017748962,385.0556014050393,361.40345006269706,290.9372420918115))
targetgene="CDK2"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
targetmat=list(c(214.03410640349637,198.60409594737004,351.0692544476919,309.2477798784222,379.47362256583193,341.7730505218324),c(212.11451800525873,237.25137948307446,233.27628091590051,321.2807674223297,361.40345006269706,279.9879910453455),c(181.40110363345656,195.38348898606134,267.92127313113326,285.18180479060726,269.24557029670933,282.33425912673107),c(222.67225419556573,346.7520161675704,206.7151202175554,276.758713509872,352.36836381112965,330.82379947536634),c(267.78258155415017,172.8392402569004,205.56028714371433,373.02261386113184,393.9297605683398,231.49845069671025),c(262.98361055855605,225.4424872916092,273.69543850033875,388.6654976682116,189.73681128291597,389.48050151000575))
targetgene="CDK4"
collabel=c("Control_Rep1","Control_Rep2","Control_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

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
Sweave("Control_vs_Condition2_Endothelial_summary.Rnw");
library(tools);

texi2dvi("Control_vs_Condition2_Endothelial_summary.tex",pdf=TRUE);

