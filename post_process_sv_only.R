#histogram of raw VAFs
args<-commandArgs(trailingOnly = TRUE)

if(length(args)==0|args[1]=="-h")
{
  print("Usage: Rscript <workingdir> <sample id> <batch no>")
}

setwd(args[1])
#setwd("/scratch/pancan/clustout_sanger/a330a96e-9897-4605-b5f1-5b5ef45cd365_sv_only/")
id<-args[2]
#id<-"a330a96e-9897-4605-b5f1-5b5ef45cd365"
batch<-args[3]

suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
convert_vafs_sv<-function(dat,p)
{
  #(1/p$purity)*(dat$adjusted_vaf)*(dat$prop_chrs_bearing_mutation*dat$most_likely_variant_copynumber)
  prob<-(dat$most_likely_variant_copynumber*dat$prop_chrs_bearing_mutation*p$purity)/(2*(1-p$purity)+dat$most_likely_variant_copynumber*p$purity)
  (1/prob)*dat$adjusted_vaf
}

#svs<-read.table(paste(id,"_filtered_adjusted_svs.tsv",sep=""),header=T,sep="\t",stringsAsFactors=F)
svs<-read.table(paste(id,"_filtered_svs.tsv",sep=""),header=T,sep="\t",stringsAsFactors=F)
#svs<-svs[order(svs$bp1_chr,svs$bp1_pos,svs$bp2_chr,svs$bp2_pos),]
mlcn<-read.table(paste(batch,"/",id,"_most_likely_copynumbers.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
pp<-read.table("purity_ploidy.txt",header=T,sep="\t",stringsAsFactors=F)
dat<-merge(svs,mlcn,by.x=c(2,3,5,6),by.y=c(1,2,3,4))
certain<-read.table(paste(batch,"/",id,"_cluster_certainty.txt",sep=""),sep="\t",header=T)
dat<-merge(dat,certain,by.x=c(1,2,3,4),by.y=c(1,2,3,4))

dat<-cbind(dat,CCF=convert_vafs_sv(dat,pp))
dat<-dat[order(dat$CCF),]
dat$CCF[dat$CCF>2]<-2

write.table(dat[,c(1,2,3,4,50)],paste0(id,"_CCF.txt"),row.names = F,sep="\t",quote = F)

pdat<-data.frame(Class=dat$classification,CCF=dat$CCF,VAF=dat$adjusted_vaf,cluster=as.character(round(dat$average_proportion*(1/pp$purity),2)))

sv_clust<-read.table(paste(batch,"/",id,"_subclonal_structure.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
sv_clust<-sv_clust[sv_clust$n_ssms>1,]#&sv_clust$n_ssms/sum(sv_clust$n_ssms)>0.01,]

plot1<-ggplot(pdat,aes(x=as.numeric(VAF),fill=Class,colour=Class))+
  geom_bar(alpha=0.3,binwidth=0.05)+xlab("Raw VAF")+xlim(0,2)
plot2<-ggplot(pdat,aes(x=as.numeric(CCF),fill=Class,colour=Class))+xlim(0,2)+
  geom_bar(alpha=0.3,binwidth=0.05)+xlab("CCF")+
  #geom_vline(xintercept=as.numeric(vals),colour="red")
  geom_vline(xintercept=1/pp$purity*as.numeric(sv_clust$proportion[sv_clust$n_ssms/(sum(sv_clust$n_ssms))>0.1&sv_clust$n_ssms>2]),colour="blue",size=1)+
  geom_vline(xintercept=1/pp$purity*as.numeric(sv_clust$proportion[sv_clust$n_ssms/(sum(sv_clust$n_ssms))<0.1|sv_clust$n_ssms<=2]),colour="red",lty=2)
plot3<-ggplot(pdat,aes(x=CCF,y=VAF,fill=cluster,colour=cluster))+geom_point(size=3)+xlim(0,2)+ylim(0,1)

bb<-read.table(paste("../",id,"_subclones.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
clon<-bb[is.na(bb$nMaj2_A),]
subclon<-bb[!is.na(bb$nMaj2_A),]

#plot3<-ggplot(bb,aes(x=as.numeric(frac1_A)))+geom_density(alpha=0.3)+xlab("Battenberg CN")

pdf(paste(id,batch,"_VAF_CCF.pdf",sep=""),height=11)
grid.arrange(plot1,plot2,plot3)
dev.off()

suppressMessages(library(circlize))
pdat<-clon[,c("chr","startpos","endpos")]
pdat<-cbind(pdat,value=clon$nMaj1_A+clon$nMin1_A)
pdat$chr<-paste("chr",pdat$chr,sep="")
colnames(pdat)<-c("chr","start","end","value")
pdat$value[pdat$value>6]<-6

pdat2<-subclon[,c("chr","startpos","endpos")]
pdat2<-cbind(pdat2,value=apply(cbind(subclon$nMaj2_A+subclon$nMin2_A,subclon$nMaj1_A+subclon$nMin1_A),1,mean))
pdat2$chr<-paste("chr",pdat2$chr,sep="")
# pdat3<-subclon[,c(2,3,4)]
# pdat3<-cbind(pdat3,value=subclon$nMaj1_A+subclon$nMin1_A)
# pdat3$chr<-paste("chr",pdat3$chr,sep="")
#pdat2<-rbind(pdat2,pdat3)
pdat2$value[pdat2$value>6]<-6
colnames(pdat2)<-c("chr","start","end","value")

pdat<-list(pdat,pdat2)
colours<-c("#0000FF80","#FF000080","darkgreen","#0000FF40","#FF000040","#00FF0040")


res_svs<-read.table(paste(batch,"/",id,"_cluster_certainty.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)

pdf(paste(id,batch,"_circos.pdf",sep=""),height=12,width=12)
par(mar = c(1, 1, 1, 1))
circos.initializeWithIdeogram(plotType = c("axis","labels"))
circos.genomicTrackPlotRegion(pdat,ylim=c(0,6),
                              panel.fun=function(region,value,...){
                                i=getI(...)
                                circos.genomicLines(region,value,type="segment",lwd=3,col=colours[i],...)
                              })

for(j in 1:nrow(dat))
{
      x<-dat[j,]
      ccf<-x$average_proportion*(1/pp$purity)
      if(ccf>1.2)
      {
        lcol=colours[1]
      }
      else if(ccf<1.2&ccf>0.9)
      {
        lcol=colours[1]
      }
      else
      {
        lcol=colours[2]
      }
      circos.link(paste("chr",as.character(x[1]),sep=""),as.numeric(x[2]),paste("chr",as.character(x[3]),sep=""),as.numeric(x[4]),col=lcol,lwd=2)
}

dev.off()
tabout<-c()
for(i in 1:nrow(sv_clust))
{
  curr<-as.numeric(sv_clust[i,c(1,3),])
  curr<-c(curr,curr[2]*1/pp$purity)
  curr<-c(curr,sum(res_svs$most_likely_assignment==curr[1]))
  tabout<-rbind(tabout,curr)
}

tabout<-tabout[order(as.numeric(tabout[,3]),decreasing=TRUE),]
tabout<-rbind(c("cluster","proportion","CCF","SVs"),tabout)


pdf(paste(id,batch,"_table.pdf",sep=""),height=2,width=3.5)
grid.table(tabout,rows=c())
dev.off()
