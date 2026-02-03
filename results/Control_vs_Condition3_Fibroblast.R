pdf(file='Control_vs_Condition3_Fibroblast.pdf',width=4.5,height=4.5);
gstable=read.table('Control_vs_Condition3_Fibroblast.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("HGF","HBEGF","TIE1","IGFBP5","PDGFB","ANGPT1","EGFR","TEK","COL5A1","IGF2")
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
targetmat=list(c(1630.5315877884625,1616.3703072626395,1622.3907828738595,444.19584102347227,536.3153046552293,443.53030540075133),c(1268.0815013448732,1709.2424398315197,1332.4578292804322,413.8407721608204,477.18576429672237,415.12688682828946),c(2075.9402771237724,1635.142334058477,2106.199875648302,648.5866380319948,357.8893232225418,650.0013096390321),c(1474.4904061587974,1590.6822705946513,1380.0456088975905,399.6750733582495,683.6204753729131,400.92517754205846),c(1978.167638001387,946.5053510743329,2037.4619717568512,667.8115149783409,492.74616965422416,668.5727756287188),c(1050.8089699617954,1636.1303354687843,1036.356089440336,246.88789341623513,657.6864664437435,245.7988145693819))
targetgene="HGF"
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
targetmat=list(c(1261.1682844372299,1359.489940582758,1328.93280856805,221.5920026973586,349.5904403652075,221.76515270037567),c(1846.8165167561629,1295.2698489127874,1825.07947383583,645.5511311457296,254.153287505863,646.7239921114403),c(1895.2090351096665,1482.0021154608553,1881.4798052339436,274.2074553926218,421.1683050097158,273.1097939659799),c(1876.4445892174915,1703.3144313696764,1879.7172948777525,519.0716775513469,580.9218000134011,520.0010477112257),c(1415.2342612361397,1916.7227359960395,1391.5019262128321,501.8704718625108,370.3376475085432,501.4295817215391),c(1763.857913864442,1626.250321365712,1735.1914456700865,644.5392955169746,740.6752950170865,645.6315529355765))
targetgene="HBEGF"
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
targetmat=list(c(1372.7673573749018,1339.7299123766131,1370.3518019385397,487.7047730599399,449.1770346532191,488.32031161117203),c(1855.7049384945615,1632.1783298275552,1853.2796395348867,640.4919530019542,386.9354132232118,642.3542354079847),c(1517.944912435413,2330.6953269147716,1517.5214166804922,338.9649356329458,502.08241286872527,339.748583693679),c(1612.7547443116653,2029.3548967710644,1609.1719552024267,543.3557326414683,729.2643310882518,545.1271487560958),c(1507.0812858662591,1314.041875708625,1505.1838441871548,509.96515689255125,454.36383643905305,511.26153430431435),c(2605.2951717661804,1545.2342057205185,2600.584030559891,547.4030751564886,668.0600700154113,548.4044662836876))
targetgene="TIE1"
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
targetmat=list(c(1616.7051539731756,1690.4704130356822,1689.3661764091194,486.6929374311849,433.6166292957173,486.1354332594442),c(1980.1428428321424,1284.401833399408,1940.5239021663435,448.2431835384925,511.4186560832264,448.99250128007094),c(1334.2508631751743,1774.4505329117974,1236.40101486802,428.0064709633912,607.8931692997377,427.14371776279256),c(1623.6183708808192,1273.5338178860284,1702.5850040805522,502.88230749126586,893.1672675206042,503.6144600732669),c(1579.1762621888258,2511.499585000996,1609.1719552024267,599.0066922229968,679.471033944246,599.7491075492918),c(1294.7467665600693,1607.4782945698744,1214.3696354156318,394.6158952144742,477.18576429672237,393.278103311011))
targetgene="IGFBP5"
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
targetmat=list(c(1963.353601770723,1351.5859293003,1961.6740264406362,351.10696317800654,384.8606925088783,350.6729754523182),c(1454.7383578512447,1649.9623552130856,1454.0710438576143,482.6455949161646,807.0663578757609,480.67323738012465),c(1413.2590564053844,1590.6822705946513,1411.7707953090294,471.5154029998589,471.9989625108884,471.93372397321326),c(1756.9446969567987,1295.2698489127874,1752.816549231997,544.3675682702234,481.3352057253895,541.8498312285042),c(1826.0768660332326,1856.454649967298,1820.6731979453523,585.852829049181,703.3303221590821,585.5473982630609),c(1059.697391700194,1553.1382170029763,1059.2687240708196,364.26082635182235,492.74616965422416,363.7822455626852))
targetgene="PDGFB"
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
targetmat=list(c(1633.4943950345953,1424.6980336630356,1632.0845898329103,688.0482275534422,297.7224225068681,553.8666621630073),c(1537.6969607429655,2004.6548615133836,1536.027775420498,509.96515689255125,572.6229171560668,633.6147220010733),c(2467.0308336133126,1556.102221233898,2461.3457124207985,394.6158952144742,392.12221500904576,372.5217589695966),c(1843.8537095100298,1697.3864229078329,1842.7045773977404,669.8351862358511,470.9616021537216,722.1022952460509),c(1307.5855979599783,1221.1697431397447,1308.663939471853,420.92362156210584,447.1023139388855,596.4717900217001),c(2050.262614323954,1145.0936345460875,2050.680799428284,641.5037886307093,391.08485465187897,572.4381281526938))
targetgene="ANGPT1"
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
targetmat=list(c(1419.1846708976502,1782.3545441942554,1399.433222815692,502.88230749126586,611.0052503712379,504.70689924913086),c(2181.6137355691785,1863.3706598394488,2112.368661894971,483.65743054491963,578.8470792990676,482.8581157318525),c(1354.0029114827269,1552.1502155926692,1278.7012634166051,386.5212101844337,601.6690071567369,386.7234682558275),c(1589.0522863426022,1464.218090075325,1570.3967273662236,616.2078979118328,560.1745928700653,617.2281343631146),c(1809.2876249718129,1147.069637366702,1788.066756355818,534.2492119826728,506.23185429739243,534.2027569974567),c(1822.126456371722,1506.7021507185361,1922.898798604433,540.3202257552032,742.75001573142,539.6649528767763))
targetgene="EGFR"
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
targetmat=list(c(1783.6099621719948,2297.103278964326,1781.0167149310537,619.243404798098,642.1260610862416,619.4130127148424),c(1341.1640800828177,1623.2863171347901,1339.5078707051964,491.75211557496016,371.37500786571,492.69006831462775),c(2544.0638220127676,1983.9068318969316,2543.3024439836818,881.3088326456591,464.7374400107209,881.5984149221831),c(1839.9032998485195,963.301375049556,1840.0608118634539,329.8584149741502,544.6141875125636,329.9166311109037),c(1416.2218636515174,1646.0103495718565,1416.1770711995068,547.4030751564886,639.0139800147413,547.3120271078237),c(1476.4656109895527,1138.1776246739369,1474.3399129538116,292.42049671021294,352.70252143670785,291.6812599556665))
targetgene="TEK"
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
targetmat=list(c(2086.8039036929263,1774.4505329117974,2080.643475483532,867.1431338430882,746.8994571600872,869.58158398768),c(1903.1098544326876,1704.3024327799835,1902.629929508236,768.9950778538472,1112.0503028827966,770.1696189840633),c(1527.8209365891892,1319.9698841704685,1524.5714581052564,485.6811018024298,740.6752950170865,486.1354332594442),c(1335.2384655905519,1586.7302649534224,1331.5765741023365,534.2492119826728,976.1560960939473,534.2027569974567),c(1706.5769737725398,1340.7179137869205,1702.5850040805522,684.000885038422,631.7524575145737,682.7744849149498),c(1496.2176592971052,1680.59039893261,1492.8462716938175,515.0243350363265,703.3303221590821,513.4464126560422))
targetgene="COL5A1"
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
targetmat=list(c(1178.2096815455093,1633.1663312378626,1142.9879659898945,379.43836078314825,369.30028715137644,380.16883320064403),c(1639.4200095268611,2091.5989856204205,1651.4722037510119,705.2494332422782,514.5307371547267,705.7157076080921),c(1368.8169477133913,1850.5266415054546,1355.3704639109158,441.16033413720703,650.424943943576,440.2529878731596),c(1811.2628298025681,1440.5060562279514,1840.9420670415493,385.5093745556786,693.994078944581,384.5385899040997),c(1416.2218636515174,2102.4670011338,1416.1770711995068,533.2373763539177,326.76851250753816,533.1103178215927),c(1664.1100699113017,1062.1015160802797,1654.1159692852984,401.6987446157596,780.0949885894244,400.92517754205846))
targetgene="IGF2"
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
targetgenelist=c("FADD","SURVIVIN","DIABLO","BIRC5","XIAP","FLIP","CFLAR","TRAF2","SCN7A","SMAD2")
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
targetmat=list(c(1332.275658344419,1266.6178080138777,1358.8954846232978,2188.600464997199,3075.7734589995257,2188.155669255431),c(1786.5727694181276,1444.4580618691803,1696.4162178338836,3097.2288596192448,2632.8205864893075,3099.2499419259398),c(1082.4122472538793,2399.855425636278,1125.362862427984,2046.9434769714906,2657.71723506131,2047.2310155689854),c(1398.44502017472,1242.905774166504,1395.9082021033098,2567.0269901515926,2703.361090776649,2568.324502456075),c(1607.816732234777,1343.6819180178422,1637.3721209014834,2753.204745842524,2570.5789650593,2750.7618448253497),c(1565.3498283735391,1353.5619321209144,1539.5527961328803,2946.465350934741,2437.7968393419515,2953.9555315360385))
targetgene="FADD"
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
targetmat=list(c(1496.2176592971052,1854.4786471466837,1560.7029204071728,3060.802776984063,2851.7036218514995,3062.1070099465664),c(1687.8125278803648,1401.9740012259692,1692.009941943406,3481.7263985461686,2655.642514346977,3483.7885318300396),c(1796.4487935719037,1233.0257600634316,1758.9853354786658,2285.736685357685,2291.5290289814343,2286.4751950831837),c(1191.0485129454185,1951.3027853567928,1246.9760770051664,2475.961783563637,2237.586290408761,2476.559611683506),c(1267.0938989294957,2361.323370634296,1182.6444490041931,2243.2395889499726,2091.318480048244,2244.962506400355),c(1317.4616221137546,1396.0459927641257,1361.5392501575843,3009.1991599175544,2339.247605411106,3008.5774903292345))
targetgene="SURVIVIN"
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
targetmat=list(c(1859.655348156072,1070.0055273627374,1905.2736950425226,3302.6314922565225,3091.3338643570273,3304.6285069883565),c(1759.9075042029317,1881.154685224979,1762.5103561910478,3958.3009796898027,2762.4906311351556,3963.3693300343),c(1628.5563829577072,1435.5660491764152,1589.7843412843251,2545.778441947736,2843.404738994165,2545.3832797629325),c(1428.0730926360488,1486.9421225123915,1328.0515533899545,2154.198053619527,2384.891461126445,2153.1976156277856),c(1305.6103931292232,1970.0748121526303,1306.901429115662,2000.3990380487576,2704.398451133816,2003.5334485344285),c(1849.7793240022957,1111.5015865956416,1783.6604804653405,2749.1574033275037,3435.737502936401,2748.5769664736217))
targetgene="DIABLO"
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
targetmat=list(c(2014.7089273703593,1162.8776599316177,1985.4679162492152,2912.0629395570686,2980.336306140181,2913.5352820290736),c(1447.8251409436014,1687.5064088047607,1411.7707953090294,2816.950390454093,2436.7594789847844,2818.4930737289123),c(1520.9077196815458,1752.7145018850383,1450.5460231452323,2670.234224284609,2901.4969189955054,2674.291102514875),c(1076.4866327616137,1162.8776599316177,1016.087220344139,3190.3177374647107,3080.9602607853594,3188.829954346781),c(1887.3082157866454,1974.0268177938592,1845.348342932027,2530.60090751641,3477.2319172230727,2531.181570476702),c(1312.5236100368666,1441.4940576382585,1378.2830985413993,2615.5951003318355,2334.0608036252725,2619.669143721679))
targetgene="BIRC5"
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
targetmat=list(c(1611.7671418962875,2124.2030321605594,1607.4094448462356,2951.524529078516,2367.2563350546097,2953.9555315360385),c(1943.6015534631701,1587.7182663637298,1928.1863296730062,3464.5251928573325,2874.525549709169,3467.401944192081),c(1880.394998879002,1441.4940576382585,1993.3992128520752,2516.4352087138395,2441.9462807706186,2515.8874220146067),c(1466.5895868357763,1682.5664017532245,1420.5833470899845,2551.8494557202666,2442.983641127785,2550.845475642252),c(1113.027922130586,1475.0861055887046,1137.7004349213214,3291.501300340217,3668.1062229417616,3293.7041152297174),c(1439.9243216205805,1983.9068318969316,1446.1397472547546,2232.109397033667,2570.5789650593,2237.3154321693073))
targetgene="XIAP"
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
targetmat=list(c(1857.6801433253167,1619.3343114935612,1852.3983843567912,3112.4063940505707,2128.6634529062485,3113.451651212171),c(1757.9322993721764,1889.058696507437,1651.4722037510119,2078.3103814628976,2158.746903264085,2080.004190844903),c(1348.077296990461,1223.1457459603591,1372.1143122947308,2082.357723977918,2381.779380054945,2084.3739475483585),c(1911.9982761710862,1250.809785448962,1823.3169634796388,2520.4825512288594,3108.968990428863,2525.7193745973823),c(1211.7881636683485,1463.2300886650178,1205.5570836346767,2981.879597941168,2329.9113621966053,2981.2665109326367),c(1323.3872366060205,1230.06175583251,1451.4272783233278,3907.7091982520496,1964.7605164738961,3910.9322495928322))
targetgene="FLIP"
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
targetmat=list(c(1711.5149858494278,2005.642862923691,1744.0039974510419,2645.9501691944874,2546.719676844464,2646.9801231182773),c(2054.2130239854646,1830.76661329931,1983.7054058930244,3725.5787850761385,2155.634822192585,3729.5873463994217),c(1618.6803588039309,1565.9822353369705,1508.7088648995368,2581.1926889541633,2429.497956484617,2581.433772566442),c(1364.8665380518808,1444.4580618691803,1234.638504511829,3242.933190159974,2647.3436314896426,3241.2670347882495),c(1771.7587331874631,1818.910596375623,1778.3729493967671,2709.6958138060563,2633.857946846474,2712.5264736701124),c(1648.3084312652597,1275.5098207066428,1671.741072847209,1778.807035351399,2174.3073086215873,1776.3060999547333))
targetgene="CFLAR"
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
targetmat=list(c(1585.1018766810917,1466.1940928959395,1513.9963959681102,2476.973619192392,2127.6260925490815,2478.7444900352334),c(1355.9781163134821,2029.3548967710644,1306.0201739375664,2181.517615595914,3106.8942697145294,2184.8783517278393),c(1567.3250332042944,1425.6860350733427,1631.2033346548146,1838.5053374479478,2889.048594709504,1840.7600113307044),c(1397.4574177593424,1475.0861055887046,1389.7394158566412,3233.8266695011785,2770.78951399249,3236.897278084794),c(1446.8375385282238,2042.1989151050586,1576.5655136128923,2401.0859470357623,2192.9797950505895,2401.1813085488952),c(1260.1806820218524,1218.2057389088232,1204.6758284565813,2771.417787160115,3037.3911257843547,2771.518189166764))
targetgene="TRAF2"
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
targetmat=list(c(1462.6391771742658,2212.1351576779034,1458.477319748092,2489.1156467374526,1676.374337181529,2488.576442618009),c(1292.771561729314,1992.7988445896967,1290.157580731847,2382.8729057181713,2468.9176500569547,2383.7022817350726),c(1448.812743358979,1469.1580971268613,1446.1397472547546,1699.8838563085042,1254.1686718146464,1700.9277968201227),c(1321.4120317752652,2422.579458073345,1317.4764912528083,1582.510923372917,1745.877481111704,1584.0368050026834),c(1081.4246448385018,2018.486881257685,1080.4188483451123,1872.9077488256198,1811.2311836132114,1875.71806495835),c(1858.6677457406943,1463.2300886650178,1857.6859154253643,2926.2286383596393,1454.3792207478364,2929.9218696670323))
targetgene="SCN7A"
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
targetmat=list(c(1194.0113201915512,1229.0737544222027,1191.4570007851482,2196.69515002724,1770.7741296837066,2199.08006101407),c(1258.205477191097,1081.8615442864243,1256.669883964217,2346.446823082989,2065.3844711190745,2348.7442281074273),c(1227.5898023143905,2021.4508854886067,1223.1821871965872,1756.5466515187877,2188.8303536219223,1756.6421947891827),c(2147.0476510309613,2112.3470152368723,2142.3313379502188,2228.0620545186466,1512.4714007491766,2228.5759187623958),c(1784.5975645873723,1444.4580618691803,1777.4916942186717,1429.7237434309027,1742.7654000402035,1433.2801987334624),c(1056.7345844540612,1352.5739307106073,1054.8624481803422,2289.7840278727053,1709.5698686108663,2290.8449517866393))
targetgene="SMAD2"
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
Sweave("Control_vs_Condition3_Fibroblast_summary.Rnw");
library(tools);

texi2dvi("Control_vs_Condition3_Fibroblast_summary.tex",pdf=TRUE);

