#histogram of raw VAFs
args<-commandArgs(trailingOnly = TRUE)

if(length(args)==0|args[1]=="-h")
{
  print("Usage: Rscript <workingdir> <sample id> <run dir> <CN dir")
}

setwd(args[1])
#setwd("/scratch/pancan/consensus_CN1/output/712ba532-fb1a-43fa-a356-b446b509ceb7/")
id<-args[2]
#id<-"712ba532-fb1a-43fa-a356-b446b509ceb7"
batch<-args[3]
#batch<-"best_run_coclus"
cndir<-args[4]
#cndir<-"../../copy-number"

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
dat<-merge(svs,mlcn,by.x=2:7,by.y=1:6)
certain<-read.table(paste(batch,"/",id,"_cluster_certainty.txt",sep=""),sep="\t",header=T)
dat<-merge(dat,certain,by.x=1:6,by.y=1:6)

dat<-cbind(dat,CCF=convert_vafs_sv(dat,pp))
dat<-dat[order(dat$CCF),]
dat$CCF[dat$CCF>2]<-2

snvs<-read.table(paste(id,"_filtered_snvs.tsv",sep=""),header=T,sep="\t",stringsAsFactors=F)
snvs_mlcn<-read.table(paste(batch,"/snvs/",id,"_most_likely_copynumbers.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
snvs_dat<-merge(snvs,snvs_mlcn,by.x=2,by.y=2)
snvs_dat<-cbind(snvs_dat,support=snvs_dat$var,depth=(snvs_dat$var+snvs_dat$ref))
snvs_dat<-cbind(snvs_dat,adjusted_vaf=(snvs_dat$support/snvs_dat$depth))
snvs_dat<-cbind(snvs_dat,CCF=convert_vafs_sv(snvs_dat,pp))
snvs_dat<-cbind(snvs_dat,classification="SNV")
snvs_dat$CCF[snvs_dat$CCF>2]<-2

certain<-read.table(paste(batch,"/snvs/",id,"_cluster_certainty.txt",sep=""),sep="\t",header=T)
snvs_dat<-merge(snvs_dat,certain,by.x=1:2,by.y=2:1)


#write.table(dat[,c(1,2,3,4,50)],paste0(id,"_CCF.txt"),row.names = F,sep="\t",quote = F)

cdat<-dat[,c("classification","CCF","adjusted_vaf","average_proportion")]
cdat<-rbind(cdat,snvs_dat[,c("classification","CCF","adjusted_vaf","average_proportion")])

pdat<-data.frame(cdat,cluster=as.character(round(cdat$average_proportion*(1/pp$purity),2)))

sv_clust<-read.table(paste(batch,"/",id,"_subclonal_structure.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
sv_clust<-sv_clust[sv_clust$n_ssms>1,]#&sv_clust$n_ssms/sum(sv_clust$n_ssms)>0.01,]

plot1<-ggplot(pdat,aes(x=as.numeric(adjusted_vaf),fill=classification,colour=classification))+
  geom_histogram(alpha=0.3,binwidth=0.05)+xlab("Raw VAF")+xlim(0,2)
plot2<-ggplot(pdat,aes(x=as.numeric(CCF),fill=classification,colour=classification))+xlim(0,2)+
  geom_histogram(alpha=0.3,binwidth=0.05)+xlab("CCF")+
  #geom_vline(xintercept=as.numeric(vals),colour="red")
  geom_vline(xintercept=1/pp$purity*as.numeric(sv_clust$proportion[sv_clust$n_ssms/(sum(sv_clust$n_ssms))>0.1&sv_clust$n_ssms>2]),colour="blue",size=1)+
  geom_vline(xintercept=1/pp$purity*as.numeric(sv_clust$proportion[sv_clust$n_ssms/(sum(sv_clust$n_ssms))<0.1|sv_clust$n_ssms<=2]),colour="red",lty=2)
plot3<-ggplot(pdat,aes(x=CCF,y=adjusted_vaf,fill=cluster,colour=cluster))+geom_point(size=3)+xlim(0,2)+ylim(0,1)+ylab("VAF")

bb<-read.table(paste(cndir,"/",id,".consensus.20161103.somatic.cna.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
bb<-bb[bb$star==3,]
#plot3<-ggplot(bb,aes(x=as.numeric(frac1_A)))+geom_density(alpha=0.3)+xlab("Battenberg CN")

pdf(paste(id,"_",batch,"_VAF_CCF.pdf",sep=""),height=11)
grid.arrange(plot1,plot2,plot3)
dev.off()

suppressMessages(library(circlize))
pdat<-bb[,c("chromosome","start","end","total_cn")]
#pdat<-cbind(pdat,value=clon$nMaj1_A+clon$nMin1_A)
pdat$chromosome<-paste("chr",pdat$chromosome,sep="")
colnames(pdat)<-c("chr","start","end","value")
pdat$value[pdat$value>6]<-6

colours<-c("#0000FF80","#FF000080","darkgreen","#0000FF40","#FF000040","#00FF0040")

res_snvs<-read.table(paste(batch,"/snvs/",id,"_cluster_certainty.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)

psnvs<-list()
count<-1
for(i in 1:nrow(sv_clust))
{
  curr<-res_snvs[res_snvs$most_likely_assignment==sv_clust[i,1],]
  if(nrow(curr)>1)
  {
    bed<-curr[,c(1,2)]
    bed[,1]<-paste("chr",bed[,1],sep="")
    bed<-cbind(bed,bed[,2]+1,as.numeric(sv_clust[i,3]))
    colnames(bed)<-c("chr","start","end","value")
    psnvs[[count]]<-bed
    count<-count+1
  }
}


res_svs<-read.table(paste(batch,"/",id,"_cluster_certainty.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)

pdf(paste(id,"_",batch,"_circos.pdf",sep=""),height=12,width=12)
par(mar = c(1, 1, 1, 1))
circos.initializeWithIdeogram(plotType = c("axis","labels"))
circos.genomicDensity(psnvs,col=colours,overlap=FALSE)
circos.genomicTrackPlotRegion(pdat,ylim=c(0,6),
                              panel.fun=function(region,value,...){
                                i=getI(...)
                                circos.genomicLines(region,value,type="segment",lwd=3,col=colours[i],...)
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(1,1), lty = 2, col = "#00000040")
                                circos.lines(cell.xlim, c(2,2), lty = 2, col = "#00000040")
                                circos.lines(cell.xlim, c(3,3), lty = 2, col = "#00000040")
                                circos.lines(cell.xlim, c(4,4), lty = 2, col = "#00000040")
                                circos.lines(cell.xlim, c(5,5), lty = 2, col = "#00000040")
                              })

for(j in 1:nrow(dat))
{
      x<-dat[j,]
      ccf<-x$average_proportion*(1/pp$purity)
      if(ccf>1.2) {
        lcol=colours[1]
      } else if(ccf<1.2 & ccf>0.9) {
        lcol=colours[1]
      } else {
        lcol=colours[2]
      }
      circos.link(paste("chr",as.character(x[1]),sep=""),
                  as.numeric(x[2]),
                  paste("chr",as.character(x[4]),sep=""),
                  as.numeric(x[5]),col=lcol,lwd=2)
}

dev.off()
tabout<-c()
for(i in 1:nrow(sv_clust))
{
  curr<-as.numeric(sv_clust[i,c(1,3),])
  curr<-c(curr,curr[2]*1/pp$purity)
  curr<-c(curr,sum(res_svs$most_likely_assignment==curr[1]))
  curr<-c(curr,sum(res_snvs$most_likely_assignment==curr[1]))
  tabout<-rbind(tabout,curr)
}

tabout<-tabout[order(as.numeric(tabout[,3]),decreasing=TRUE),]
tabout<-rbind(c("cluster","proportion","CCF","SVs","SNVs"),tabout)


pdf(paste(id,"_",batch,"_table.pdf",sep=""),height=2,width=3.5)
grid.table(tabout,rows=c())
dev.off()
