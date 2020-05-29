R版本：3.6.2
##原始数据的质量分析
setwd("")
library(affyPLM)
Data<-ReadAffy()
Pset<-fitPLM (Data) #对数据集进行回归计算
image(Data[,1]) #灰度图
image(Pset,type="weights",which=1,main="Weights")#权重图
image(Pset,type="resids",which=1,main="Residuals")#残差图

#质量控制:相对对数表达（RLE）
library(affyPLM)
library(RColorBrewer) 
colors<-brewer.pal(12,"Set3")
Mbox(Pset,col=colors,main="RLE",las=3) #绘制RLE图
boxplot(Pset,col=colors,main="NUSE",las=3) #绘制NUSE图

#RMA法预处理baseline样本
setwd("")
library(affyPLM)
library(affy)
Data<-ReadAffy()
eset.rma<-rma(Data) #调用rma函数用RMA算法预处理数据
baseline_exprs<-exprs(eset.rma)
probeid<-rownames(baseline_exprs)
baseline_exprs<-cbind(probeid,baseline_exprs)
write.table(baseline_exprs,file="baseline.expres.txt",sep='\t',quote=F,row.names=F)

#RMA法预处理treated样本
setwd("")
Data<-ReadAffy()
eset.rma<-rma(Data)
treated_exprs<-exprs(eset.rma)
probeid<-rownames(treated_exprs)
treated_exprs<-cbind(probeid,treated_exprs)
write.table(treated_exprs,file="treated.expres.txt",sep='\t',quote=F,row.names=F)

#合并以上两个数据
setwd(" ")
baseline_exprs<-read.table("baseline.expres.txt",header=T,sep="\t")
treated_exprs<-read.table("treated.expres.txt",header=T,sep="\t")
probe_exprs<-merge(baseline_exprs,treated_exprs,by="probeid")
write.table(probe_exprs,file="skin.probeid.exprs.txt",sep='\t',quote=F,row.names=F)             ###预处理结束

#Probe ID转换为Gene symbol
setwd("")
probe_exp<-read.table("skin.probeid.exprs.txt",header=T,sep="\t",row.names=1)
probeid_geneid<-read.table("GPL571-17391.txt",header=T,sep="\t")
probe_name<-rownames(probe_exp)
loc<-match(probeid_geneid[,1],probe_name)
probe_exp<-probe_exp[loc,]
raw_geneid<-as.numeric(as.matrix(probeid_geneid[,3]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
rownames(gene_exp_matrix)<-levels(geneidfactor)
geneid<-rownames(gene_exp_matrix)
gene_exp_matrix2<-cbind(geneid,gene_exp_matrix)
write.table(gene_exp_matrix2,file="Skin.Isotretinoin_Treatment.7weeks.geneid.exprs.txt",sep='\t',quote=F,row.names=F)
loc<-match(rownames(gene_exp_matrix),probeid_geneid[,3])
rownames(gene_exp_matrix)=probeid_geneid[loc,2]
genesymbol<-rownames(gene_exp_matrix)
gene_exp_matrix3<-cbind(genesymbol,gene_exp_matrix)
write.table(gene_exp_matrix3,file="Skin.Isotretinoin_Treatment.7weeks.genesyb.exprs.txt",sep='\t',quote=F,row.names=F)
#补充缺失值(需要对genesyb这个文件进行处理)
library(impute)
#读取表达值,第一行会报错 需要删除没有基因名字的表达量
gene_exp_matrix<-read.table("Skin.Isotretinoin_Treatment.7weeks.genesyb.exprs.txt",header=T,sep="\t",row.names=1)
gene_exp_matrix<-as.matrix(gene_exp_matrix)
#利用KNN法补充缺失值
imputed_gene_exp<-impute.knn(gene_exp_matrix,k=10,rowmax=0.5,colmax=0.8,maxp=3000,rng.seed=362436069)
GeneExp<-imputed_gene_exp$data
genesymbol<-rownames(GeneExp)
GeneExp<-cbind(genesymbol,GeneExp)
write.table(GeneExp,file="gene.Skin.Isotretinoin_Treatment.7weeks.gene.exprs.txt",sep='\t',quote=F,row.names=F)
## 经过预处理 基因注释 得到完整的基因表达矩阵
##做差异基因
library(limma)
rt<-read.table("gene.Skin.Isotretinoin_Treatment.7weeks.gene.exprs.txt",header=T,sep="\t",row.names="genesymbol")
#differential
class<-c(rep("baseline",8),rep("treated",8))
design<-model.matrix(~factor(class))
colnames(design)<-c("baseline","treated")
#计算均值和方差，并通过贝叶斯检验得到差异结果
fit<-lmFit(rt,design)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table 对总表进行过滤，两组之间的差异是两倍以上的
diffLab<-allDiff[with(allDiff, ((logFC>1 |logFC<(-1)) & adj.P.Val<0.05)),]
write.table(diffLab,file="diffEXp.xls",sep="\t",quote=F)

#heatmap 热图
hmExp=log10(diffExpLevel+0.00001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",height=120,width=90)
par(oma=c(3,3,3,5))
heatmap.2(hmMat,col='greenred',trace="none",cexCol=1)
dev.off()
#volcano 火山图
pdf(file="vol.pdf")
#确定x轴的最大值 x轴表示调整后的p值
xMax=max(-log10(allDiff$adj.P.Val))
yMax=max(abs(allDiff$logFC))#y轴代表logFC
#对所有基因进行打点
plot(-log10(allDiff$adj.P.Val),allDiff$logFC,xlab="adj.P.Val",ylab="logFC",main="Volcano",xlim=c(0,xMax),ylim=c(-yMax,yMax),pch=20,cex=0.4)
#将有差异的基因打点成红色
diffSub=subset(allDiff,allDiff$adj.P.Val<0.05 & abs(allDiff$logFC)>1)
points(-log10(diffSub$adj.P.Val),diffSub$logFC,pch=20,col="red",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()
#GO富集分析，可视化富集分析结果
setwd("G:\\BYSJ\\plus\\GES11792")
library（clusterprofiler）#基因富集分析用
library（org.Hs.eg.db）
gnames<-read.table('geneNames.txt',header = T,sep = '\t')
#选取基因列的所有行
b<-gnames[,1]
#利用bitr函数将基因名称转换为ENTREZID号，物种是人org.Hs.eh.db
eg<-bitr(b,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
#可能会有部分基因对应不到ENTREZID
#转换后的基因名称保存为文档
write.table(eg,file = "go_id.txt")
gene<-eg[,2]
#进行GO分析
ego_CC<-enrichGO(gene = gene,OrgDb = org.Hs.eg.db,ont = "CC",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.01,qvalueCutoff = 0.01,readable = TRUE)
write.csv(as.data.frame(ego_CC),row.names = F,file = "ego_CC.csv")
#绘图
barplot(ego_CC,drop=TRUE,title="enrichment_CC",showCategory=12)
#BP
ego_BP<-enrichGO(gene = gene,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.01,qvalueCutoff = 0.01,readable = TRUE)
write.csv(as.data.frame(ego_BP),row.names = F,file = "ego_BP.csv")
barplot(ego_BP,drop=TRUE,title="enrichment_BP",showCategory=12)
#MF
ego_MF<-enrichGO(gene = gene,OrgDb = org.Hs.eg.db,ont = "MF",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.01,qvalueCutoff = 0.01,readable = TRUE)
write.csv(as.data.frame(ego_MF),row.names = F,file = "ego_MF.csv")
barplot(ego_MF,drop=TRUE,title="enrichment_MF",showCategory=12)
##多组数据联合分析
setwd("")
options(stringsAsFactors = F)
gn10433 = read.table("gn10433.txt",header=F, sep="\t")
gn10433.Gene = gn10433$V1
gn11792 = read.table("gn11792.txt", header=F, sep="\t")
gn11792.Gene = gn11792$V1
noacne.alesion = read.table("noacne.alesion.txt", header=F, sep="\t")
noacne.alesion.Gene = noacne.alesion$V1
normal.alesion = read.table("normal.alesion.txt", header=F, sep="\t")
normal.alesion.Gene = normal.alesion$V1

####交集
intersect(intersect(intersect(gn10433.Gene, gn11792.Gene), noacne.alesion.Gene), normal.alesion.Gene)
intersect(noacne.alesion.Gene, gn11792.Gene)
####频数统计
all.genes = c(gn10433.Gene,gn11792.Gene,noacne.alesion.Gene,normal.alesion.Gene)
tbl = as.data.frame(table(all.genes))
tbl[tbl$Freq == 4,]
tbl[tbl$Freq == 3,]
tbl[tbl$Freq == 2,]

####可视化，韦恩图
install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(list("one week" = gn10433.Gene,
                  "eight weeks" = gn11792.Gene,
                  "noacne-alesion" = noacne.alesion.Gene,
                  "normal-alesion" = normal.alesion.Gene) ,
             height=5000,
             width=5200,
             resolution=500,
             imagetype="tiff",
             filename="VennPlot.tiff",
             col="transparent",
             fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),alpha = 0.50,
             label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
             cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
             cex=1.5,
             cat.cex=1.4
)
