pdf(file='Control_vs_Condition1_Fibroblast.pdf',width=4.5,height=4.5);
gstable=read.table('Control_vs_Condition1_Fibroblast.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("CNN1","MMP2","FN1","COL3A1","TAGLN","MYH11","VIM","PLAU","SMAD3","COL1A1")
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
targetmat=list(c(2650.8043228174333,1808.4615061274571,2651.188423344609,455.4209513266546,593.5587610932089,521.974774340197),c(1593.9476973804174,1602.4951610975074,1600.3600034354397,564.2776665218062,312.130038161084,441.326311352012),c(1337.915035844674,1535.122992162477,1341.959572310234,355.4504985964133,440.05218494841347,483.8907779291096),c(1561.221718086375,1330.1191066887418,1565.0452778483284,607.598196038244,444.145693645608,375.2393764033605),c(2380.3337292402007,1250.2349635229202,2385.035979285647,307.6868378475203,576.1613491301321,519.7345392571918),c(1826.879667649778,1879.6835132873464,1831.19772190729,604.2658476139027,495.3145523605398,574.6202987908177))
targetgene="CNN1"
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
targetmat=list(c(1892.3316262378628,1692.0038998254763,1899.2431687702608,473.1934762564752,293.70924902370854,318.11338178672946),c(1386.0414759829716,1711.2530909497707,1385.0263108311017,532.0649650865062,666.218540468412,470.4493674310788),c(1579.5097653389282,1616.932054440728,1652.9014244308983,373.223023526234,446.1924479942053,855.7698017079624),c(1908.694615884884,1810.3864252398866,1910.4405207856864,537.6188791270752,563.8808230385484,427.8849008539812),c(2229.216707205946,1229.0608532861963,2334.2172278310236,358.7828470207547,573.0912176072361,730.3166370596747),c(1302.3014701423338,1731.4647416302796,1190.3646527167803,439.86999201306145,335.6677131699526,356.1973781978168))
targetgene="MMP2"
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
targetmat=list(c(1506.3575763287158,1535.122992162477,1509.0585177712005,593.1580195327647,565.9275773871457,440.20619381050943),c(1885.593924618501,1230.9857723986258,1886.3231472140005,538.7296619351889,647.7977513310366,477.1700726800942),c(1391.8166487995675,940.3229864217806,1393.6396585352752,485.4120871457269,502.47819258063026,560.0587707512843),c(2457.3360334614767,1589.983186866716,2469.446786786548,406.54650776964775,392.9768349306762,733.6769896841824),c(1716.1888553316933,1908.557299973788,1717.5015322121997,408.7680733858753,547.5067882497702,424.5245482294735),c(1669.9874727989277,1685.2666829319733,1673.5734589209146,420.986684275127,565.9275773871457,460.3683095575557))
targetgene="FN1"
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
targetmat=list(c(1645.9242527297788,1449.464091659367,1649.4560853492287,380.9985031830305,566.9509545614443,485.0108954706122),c(1768.1654106810547,1265.6343164223556,1770.042953207658,632.0354178167474,717.3873991833437,537.6564199212329),c(1665.1748287850978,1443.6893343220786,1671.85078938008,474.304259064589,395.0235892792735,585.8214742058434),c(1715.2263265289275,1015.3948318065287,1721.8082060642864,516.514005772913,754.2289774580947,529.815597130715),c(1181.9853697965898,1384.0168418367662,1184.3353093238588,497.63069803497865,365.34565122461305,432.36537101999147),c(1847.092772507863,1438.877036541005,1851.0084216268892,945.2761697048367,556.717182818458,412.2032552729452))
targetgene="COL3A1"
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
targetmat=list(c(1653.6244831519064,1661.205194026605,1658.9307678238197,568.7207977542613,432.888544728323,612.704295201905),c(1496.732288301056,1548.597425949483,1501.3065048374442,664.2481192520473,522.945736066603,574.6202987908177),c(1594.9102261831833,1653.5055175768875,1593.4693252721008,448.7562544779718,650.8678828539324,424.5245482294735),c(1590.0975821693537,1625.5941904466606,1594.3306600425183,405.43572496153394,496.33792953483845,469.32924988957626),c(1554.4840164670134,1650.6181389082433,1555.5705953737374,788.655793760792,573.0912176072361,532.0558322137201),c(1447.6433193599926,2624.627209797539,1449.626418612403,393.21711407228224,472.8002545259698,287.8702081661601))
targetgene="TAGLN"
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
targetmat=list(c(1387.0040047857376,1416.7404667480666,1391.9169889944405,480.96895591327177,369.4391599218076,418.9239605219607),c(1425.5051568963756,1730.502282074065,1429.8157188928042,367.66910948566505,655.9847687254256,519.7345392571918),c(1417.804926474248,1897.0077852992115,1420.3410364182132,544.2835759757579,794.1406872557415,501.81265859315073),c(1524.6456235812689,1791.1372341155923,1527.1465479499648,582.0501914516268,608.9094187076884,426.76478331247864),c(1835.5424268746715,1638.106164677452,1837.2270653002115,737.5597845875576,412.4210012423503,420.04407806346325),c(1528.4957387923325,2085.6498583172965,1534.898560883721,402.10337653719256,546.4834110754716,521.974774340197))
targetgene="MYH11"
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
targetmat=list(c(1274.3881348621212,1713.1780100622,1277.3594645289327,573.1639289867164,452.3327110399971,459.2481920160531),c(2105.050491649138,1355.1430551503245,2110.270187522512,399.881810920965,483.03402626895615,543.2570076287458),c(1640.1490799131832,1373.429786718404,1643.4267419563073,549.8374900163268,692.8263470001765,459.2481920160531),c(1356.2030830972271,1487.9624739079557,1360.0476024889986,526.5110510459372,505.54832410352617,815.4455702138699),c(1838.4300132829694,1779.5877194410157,1840.6724043818808,544.2835759757579,627.3302078450638,436.84584118600173),c(1772.0155258921186,1576.50875307971,1776.9336313709969,479.85817310515796,536.2496393324852,457.007956933048))
targetgene="VIM"
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
targetmat=list(c(1308.0766429589296,1140.5145741144422,1308.3675162639574,747.5568298605817,505.54832410352617,421.1641956049658),c(1514.0578067508434,1468.7132827836615,1517.671865475374,566.4992321380337,343.8547305643417,537.6564199212329),c(1498.657345906588,1769.9631238788684,1498.722500526192,345.45345332338917,523.9691132409016,405.48255002392983),c(1763.352766667225,2144.359891246394,1768.3202836668233,358.7828470207547,364.3222740503144,689.9924055655822),c(1604.5355142108428,1184.7877137003193,1613.2800249917,878.6292012180091,504.5249469292275,445.8067815180223),c(1663.249771179566,1677.5670064822555,1668.4054502984104,572.0531461786027,532.1561306352907,362.9180834468322))
targetgene="PLAU"
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
targetmat=list(c(1449.5683769655245,1420.5903049729254,1450.4877533828205,730.8950877388748,583.3249893502225,287.8702081661601),c(1457.268607387652,1510.0990437008943,1460.8237706278287,390.9955484560546,555.6938056441593,395.40149215040674),c(1310.9642293672273,1635.2187860088077,1316.1195291977135,645.364811514113,465.63661430587933,542.1368900872432),c(1834.5798980719055,2679.487404501778,1833.7817262185422,488.7444355700683,605.8392871847925,486.13101301211475),c(1969.3339304591389,1737.2394989675681,1969.8726199444836,732.0058705469887,387.859949059183,579.100768956828),c(2044.4111770748832,1741.089337192427,2048.2540840524625,506.51696049988897,571.0444632586389,639.5871161979667))
targetgene="SMAD3"
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
targetmat=list(c(2150.289345379138,1412.8906285232076,2153.3369260433797,774.2156172553127,657.0081458997242,501.81265859315073),c(1511.1702203425455,1611.1572971034398,1513.3651916232873,447.645471669858,658.031523074023,613.8244127434076),c(1659.3996559685022,1524.5359370441151,1663.2374416759064,604.2658476139027,461.5431056086848,935.2981471546448),c(1349.4653814778655,1666.0174918076789,1352.2955895552425,420.986684275127,621.189944799272,561.1788882927868),c(1956.8210560231814,1938.3935462164443,1958.675267929058,771.9940516390851,447.2158251685039,879.2922700795164),c(1529.4582675950985,2340.701640714197,1530.5918870316343,548.726707208213,669.2886719913079,654.1486442375001))
targetgene="COL1A1"
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
targetgenelist=c("RICTOR","RHEB","PRAS40","DEPTOR","TSC1","INPP4B","DDIT4","AKT1S1","RAPTOR","INPP5D")
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
targetmat=list(c(1213.7488202878662,1501.4369076949617,1208.4526828955447,2451.4976575071382,2682.2715738367247,1610.7290246806936),c(927.8777658663782,1322.419430239024,993.9803250616241,2469.2701824369587,2516.484471600346,2608.7537541594825),c(1198.348359443611,1441.7644152096493,1242.9060737122388,2474.8240964775277,2173.6531182103026,1971.4068730445208),c(1773.9405834976506,1722.8026056243473,1799.328335401848,2915.804871298703,3910.324182995088,1869.476176767787),c(1203.1610034574408,1474.4880401209498,1186.0579788646935,2998.0027990991234,3368.9576577911093,3508.208139986045),c(1988.5845065144579,1167.4634416884544,1928.5285509644507,2523.6985400345343,2623.9390749017025,2511.303528048759))
targetgene="RICTOR"
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
targetmat=list(c(1491.9196442872264,1374.3922462746189,1478.0504660361757,2389.2938202527657,2672.0378020937383,3519.4093154010707),c(1258.0251452150999,1803.6492083463836,1304.0608424118707,3055.763505121041,2578.9104792325625,2902.224550033155),c(1683.462876037651,1834.4479141452546,1649.4560853492287,3567.8343796614986,3159.165337059889,2940.3085464442424),c(2048.261292285947,1268.5216950909999,2024.9980452511943,3013.5537584127164,2643.3832412133766,2392.5710686494867),c(1557.3716028753113,1830.5980759203958,1540.066569506225,2082.717765213359,2707.8560031941906,2491.1414123017125),c(1494.8072306955244,2015.3903107136218,1589.162651420014,3163.5094375080785,2644.4066183876753,2438.495887851092))
targetgene="RHEB"
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
targetmat=list(c(1451.4934345710565,1374.3922462746189,1480.6344703474279,3099.0840346374785,3868.3657188488437,3163.211937203254),c(1281.1258364814828,1348.4058382568214,1273.9141254472634,2161.5833445894386,2410.0532454732875,2804.774323922432),c(2076.1746275661594,1505.2867459198208,2178.315634385483,2138.2569056190487,2442.801315050844,2052.055336032706),c(1416.842397671482,1678.5294660384702,1436.7063970561428,2895.810780752655,2653.617012956363,2619.954929574508),c(1447.6433193599926,1772.8505025475126,1553.8479258329028,2920.248002531158,2719.1131521114758,2196.550498886537),c(1790.3035731446716,1806.5365870150279,1733.0055580797118,2178.245086711145,3798.776070996537,2168.5475603489726))
targetgene="PRAS40"
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
targetmat=list(c(1197.385830640845,1292.5831839963678,1161.0792705225904,3062.4282019697234,2853.175561944597,1814.5904172341611),c(1649.7743679408427,1651.580598464458,1715.7788626713648,2963.568532047596,2619.845566204508,2390.3308335664815),c(1222.4115795127598,1833.48545458904,1205.0073438138752,2684.7620472110343,2051.871234468765,2513.543763131764),c(1330.2148054225463,1749.7514731983595,1272.1914559064287,2776.957020284479,2443.824692225143,2752.128799471811),c(1928.9077207429689,1526.4608561565444,1914.7471946377732,2865.8196449335824,2573.7935933610693,2129.343446396383),c(1475.5566546402054,1538.0103708311212,1411.7276887140397,2172.691172670576,1918.8322018099425,4072.7473809033395))
targetgene="DEPTOR"
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
targetmat=list(c(1520.795508370205,1661.205194026605,1502.1678396078617,1886.109208177218,2828.61450976143,2802.5340888394267),c(1706.5635673040338,1639.0686242336667,1740.7575710134681,4073.2405573532737,2045.7309714229732,2290.6403723727526),c(1303.2639989450997,2066.400667193002,1271.3301211360113,3075.757595667089,2246.312897585506,3342.4307438436645),c(1823.99208124148,1581.3210508607835,1785.5469790751704,2798.061893638641,3029.1964359239623,2560.588699874872),c(1610.3106870274387,1519.7236392630416,1649.4560853492287,2361.5242500499207,2189.003775824782,2086.7789798192853),c(1534.2709116089284,2113.5611854475233,1522.839874097878,3291.2494604411645,2749.814467340435,2798.053618673416))
targetgene="TSC1"
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
targetmat=list(c(1118.458468814037,1772.8505025475126,1181.7513050126067,2168.248041438121,1983.3049637907566,3384.9952104207623),c(1764.315295469991,1223.286095948908,1716.6401974417822,2605.8964678349553,2662.8274075250506,2019.571927329131),c(1633.4113782938214,1655.4304366893168,1535.7598956541383,2202.6823084896487,2649.5235042591685,3065.7617110925303),c(1592.9851685776514,1383.0543822805514,1591.7466557312662,1879.4445113285353,2782.562536917991,2434.0154176850815),c(1803.778976383395,1350.330757369251,1844.1177434635504,2631.4444724215723,2142.9518029813435,2273.838609250214),c(1262.8377892289298,1519.7236392630416,1310.9515205752095,2429.2820013448622,2130.67127688976,2623.3152821990157))
targetgene="INPP4B"
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
targetmat=list(c(1082.8449031116966,2012.5029320449778,1082.6978064146113,1908.324864339494,3078.318540290297,2689.402217147667),c(1354.2780254916952,1654.467977133102,1330.7622202948087,3174.6172655892165,2035.4971996799868,1744.0230121194993),c(1586.2474669582898,1204.0369048246137,1562.4612735370763,3695.5744025945846,2993.37823482351,1923.2418187599103),c(2178.202680659351,1227.1359341737668,2263.5877766568005,2327.089982998393,1980.2348322678606,2367.92848273643),c(1917.3573751097774,1284.88350754665,1844.9790782339676,2439.2790466178863,3326.9991936448655,2826.0565572109804),c(1360.053198308291,1151.101629232804,1292.0021556260278,2673.6542191298963,2275.9908356401666,3317.788157930608))
targetgene="DDIT4"
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
targetmat=list(c(1509.2451627370135,1605.3825397661515,1472.8824574136715,2644.773866118938,1925.995842030033,2481.0603544281894),c(1411.0672248548865,2017.3152298260513,1437.5677318265602,2781.4001515169343,2698.645608625503,2494.5017649262204),c(1753.7274786395656,1249.2725039667055,1724.3922103755383,2023.8462763833284,2148.0686888528367,2912.3056079066782),c(1283.0508940870147,1461.9760658901585,1285.9728122331062,2473.7133136694138,2368.0947813270436,3886.807869013913),c(1748.9148346257357,2052.926233405996,1755.400262110563,2711.4208346057653,2135.7881627612533,2050.935218491203),c(1150.2219193053133,1549.5598855056978,1177.44463116052,2059.3913262429696,2135.7881627612533,2212.232144467573))
targetgene="AKT1S1"
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
targetmat=list(c(2099.2753188325423,1801.724289233954,2053.4220926749667,2365.967381282376,1860.49970287492,3013.1161866419097),c(1607.4231006191408,1464.8634445588025,1626.2000465479603,3939.9466203796187,2972.9106913375376,2030.7731027441569),c(1303.2639989450997,1196.3372283748959,1327.3168812131391,2712.5316174138793,1612.8424266946502,3171.052759993772),c(1458.231136190418,2124.1482405658853,1460.8237706278287,2517.0338431858518,2210.4946964850537,2231.274142673117),c(1492.8821730899924,1216.548879055405,1540.066569506225,3574.4990765101816,2391.632456335912,3443.241322578896),c(1369.6784863359505,1971.1171711277448,1416.0343625661264,2413.731042031269,1977.1647007449646,1924.361936301413))
targetgene="RAPTOR"
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
targetmat=list(c(1926.020134334671,892.2000086110446,1896.6591644590087,2099.379507335066,3336.209588213553,2091.2594499852958),c(1775.8656411031825,1282.9585884342207,1895.7978296885915,2919.1372197230444,3181.679634894459,2813.7352642544524),c(2100.2378476353083,2029.8272040568427,1981.9313067303267,3706.6822306757226,2293.388247603243,2822.696204586473),c(1583.3598805499919,1145.3268718955158,1552.9865910624853,2465.9378340126173,3012.822401135184,2867.5009062465756),c(1361.9782559138227,1344.5560000319626,1382.4423065198498,3011.332192796489,3439.5706828177154,3107.2060601281255),c(1841.317599691267,921.0737952974862,1832.0590566777073,2766.959975011455,2515.461094426047,2759.969622262329))
targetgene="INPP5D"
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
Sweave("Control_vs_Condition1_Fibroblast_summary.Rnw");
library(tools);

texi2dvi("Control_vs_Condition1_Fibroblast_summary.tex",pdf=TRUE);

