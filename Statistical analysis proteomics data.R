TMPomics.tool<-function(data=Dataset,names_col=c("SAMPLE_ID","Treatment","ind1","ind2","ind3","ind4","ind5"),groups=2,label.volcano="FALSE", p_threshold= 0.05,n4group=5,PCA_dimensions=5, fold_threshold=1, p_col="red", fold_col="blue",both_col="green", cex_volcanolab= 0.6, x_lim_volcano = c(-5,3), y_lim_volcano = c(0,3),pch_volcano = 20)
{
  library(stats)
  library(RcmdrMisc)
  library(base)
  library(car)
  library(graphics)
  dataset1 <- data
  colnames(dataset1) <- names_col
  num<-length(dataset1$SAMPLE_ID)/groups
  Molecules<-data.frame(SAMPLE=c(dataset1$SAMPLE_ID[1:num]))
  ind1<-names_col[3]
  ind2<-names_col[4]
  ind3<-names_col[5]
  ind4<-names_col[6]
  ind5<-names_col[7]
  ind6<-names_col[8]
  ind7<-names_col[9]
  ind8<-names_col[10]
  ind9<-names_col[11]
  ind10<-names_col[12]
  ind11<-names_col[13]
  ind12<-names_col[14]
  ind13<-names_col[15]
  ind14<-names_col[16]
  ind15<-names_col[17]
  ind16<-names_col[18]
  ind17<-names_col[19]
  ind18<-names_col[20]
  ind19<-names_col[21]
  ind20<-names_col[22]
  if (groups==2)
  {
    #Significance hunter cycle
    i<- 1
    while(i <= num)
    {
      i1<- c(dataset1$ind1[i],dataset1$ind2[i],dataset1$ind3[i],dataset1$ind4[i],dataset1$ind5[i],dataset1$ind6[i],dataset1$ind7[i],dataset1$ind8[i],dataset1$ind9[i],dataset1$ind10[i],dataset1$ind11[i],dataset1$ind12[i],dataset1$ind13[i],dataset1$ind14[i],dataset1$ind15[i],dataset1$ind16[i],dataset1$ind17[i],dataset1$ind18[i],dataset1$ind19[i],dataset1$ind20[i])
      i2<- i+num
      i3<- c(dataset1$ind1[i2],dataset1$ind2[i2],dataset1$ind3[i2],dataset1$ind4[i2],dataset1$ind5[i2],dataset1$ind6[i2],dataset1$ind7[i2],dataset1$ind8[i2],dataset1$ind9[i2],dataset1$ind10[i2],dataset1$ind11[i2],dataset1$ind12[i2],dataset1$ind13[i2],dataset1$ind14[i2],dataset1$ind15[i2],dataset1$ind16[i2],dataset1$ind17[i2],dataset1$ind18[i2],dataset1$ind19[i2],dataset1$ind20[i2])
      e<- normalityTest(i1,test="shapiro")$p.value  
      f<-normalityTest(i3,test="shapiro")$p.value
      if(e>=0.05 & f>=0.05)
      {
        a<- mean(i1)
        a1<-sd(log2(i1), na.rm = T)
        b<- mean(i3)
        b1<-sd(log2(i3), na.rm = T)
        if((sum(is.na(i1)) < (length(i1)-2)) & (sum(is.na(i3)) < (length(i3)-2)))
        {
          c<- t.test(i3,i1,alternative="two.sided",conf.level = .95)$p.value  
          d<- log2(b/a)
          d0<-sqrt((a1^2)+(b1^2))
          Molecules[i,"p.value"] <- c
          Molecules[i,"LOG2_FOLD_CHANGE"] <- d
          Molecules[i,"±"]<-d0
          i<-i+1
        }
        else 
        {
          i<-i+1
        }
      }
      else  
      {
        a<- median(i1)
        a1<-IQR(log2(i1),na.rm=T)
        b<- median(i3)
        b1<-IQR(log2(i3),na.rm=T)
        if((sum(is.na(i1)) < (length(i1)-2)) & (sum(is.na(i3)) < (length(i3)-2)))
        {
          c<- wilcox.test(i3,i1,alternative="two.sided",conf.level = .95)$p.value  
          d<- log2(b/a)
          d0<-sqrt((a1^2)+(b1^2))
          Molecules[i,"p.value"] <- c
          Molecules[i,"LOG2_FOLD_CHANGE"] <- d
          Molecules[i,"±"]<-d0
          i<-i+1
        }
        else 
        {
          i<-i+1
        }
      }
    }
    #Save entire molecules table
    dir.create(file.path("~","Desktop","Omics"))
    write.csv(Molecules,"~/Desktop/Omics/Molecules.csv", row.names = TRUE) #Select the path
    #save significative molecules
    interesting_molecules<-data.frame()
    i<-1
    j<-0
    while(i <= num)
    {
      if ((sum(is.na(Molecules[i,])) < 1 & Molecules$`p.value`[i]< p_threshold & Molecules$LOG2_FOLD_CHANGE[i] > fold_threshold)|(sum(is.na(Molecules[i,])) < 1 & Molecules$`p.value`[i]< p_threshold & Molecules$LOG2_FOLD_CHANGE[i] < -fold_threshold))
      {
        j<-j+1
        id<-Molecules$SAMPLE[i]
        pval<-Molecules$`p.value`[i]
        fol<-Molecules$LOG2_FOLD_CHANGE[i]  
        conf<-Molecules$`±`[i]
        interesting_molecules[j,"ID"] <- id  
        interesting_molecules[j,"p-value"] <- pval
        interesting_molecules[j,"Log2_Fold_Change"] <- fol
        interesting_molecules[j,"±"]<-conf
        i<-i+1  
      }
      else
      {
        i<-i+1
      }
    }
    write.csv(interesting_molecules,"~/Desktop/Omics/Significative Molecules.csv", row.names = TRUE) #Select the path
    #Frequency distribution of Log2 FC and p-value
    pdf(file = "~/Desktop/Omics/p-value Distribution.pdf")
    par(mar=c(5, 5, 5, 5))
    par(pin=c(3.14,3.14)) 
    hist(Molecules$`p.value`,breaks="Sturges",freq=TRUE,main="p-value Distribution",xlab = "p-value")
    dev.off()
    pdf(file = "~/Desktop/Omics/Log2 Fold Change Distribution.pdf")
    par(mar=c(5, 5, 5, 5))
    par(pin=c(3.14,3.14))
    hist(Molecules$LOG2_FOLD_CHANGE,breaks="Sturges",freq=TRUE,main="Log2 Fold Change Distribution",xlab = "Log2 Fold Change")
    dev.off()
    #PCA analysis
    library(FactoMineR)
    bordind<-n4group+2
    dataset_transposed_old<- data.frame(t(dataset1[1:num,3:bordind])) #select only individuals representing variables
    othergroup<-num+1
    dataset_transposed_young<- data.frame(t(dataset1[othergroup:length(dataset1$SAMPLE_ID),3:bordind]))
    dataset_transposed<-rbind(dataset_transposed_old,dataset_transposed_young)
    tfgh<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
    ghdf<-n4group*2
    rownames(dataset_transposed)<- tfgh[1:ghdf]
    ID<-c(dataset1$SAMPLE_ID[1:num])
    colnames(dataset_transposed)<- ID
    PCA_dataset<-PCA(dataset_transposed,axes=c(1,2),graph=F,ncp=PCA_dimensions)
    PCA_dataset1<- PCA(dataset1[1:num,3:bordind],axes=c(1,2),graph=F,ncp=PCA_dimensions)
    #Save PCA dimensional results
    summary(PCA_dataset)
    summary(PCA_dataset1)
    summary_PCA<-PCA_dataset$eig
    write.csv(summary_PCA,"~/Desktop/Omics/summary_PCA_individuals.csv", row.names = TRUE) #Select the path
    #Save PCA graphs
    pdf(file = "~/Desktop/Omics/PCA_barplot_individuals.pdf")
    par(mar=c(5, 5, 5, 5))
    par(pin=c(3.14,3.14)) 
    barplot(summary_PCA[,2],ylim=c(0,100), cex.axis = 0.8, las=2)
    dev.off()
    pdf(file = "~/Desktop/Omics/PCA_individuals.pdf")
    par(mar=c(5, 5, 5, 5))
    par(pin=c(3.14,3.14)) 
    plot(PCA_dataset,choix = "ind",title="PCA Analysis of individuals",graph.type = "classic", cex.lab=1.4,cex=0.8,lwd=0.4)
    dev.off()
  }
}

    