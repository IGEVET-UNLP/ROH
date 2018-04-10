library(GenomicRanges)
library(BSgenome.Btaurus.UCSC.bosTau8)
library(stringr)
library(data.table)
library(gtools)
library(reshape2)
library(ggplot2)
library(plyr)


workingPath <- "Escritorio/Estu/pedmap_florida/ROHs_cabras1Megagap"

##Importar salidas del cgaTOH, filtrar datos y generar lista de listas por cromosoma y muestra

ff <- list.files(
		path = workingPath ,
		full.names = TRUE, 
		pattern = ".homozygousruns")

myfilelist <- lapply(
	ff, 
	read.table, 
	fill = TRUE, 
	stringsAsFactors = F, 
	sep = ",", 
	header = T)

names(myfilelist) = unlist(lapply(ff, function(x) str_match(x,"[0-9]{1,2}_[0-9]{2}_[0-9]{7,8}")))
vector = rep(1:29, each=5)
crom = split(myfilelist,vector)
rm(vector)

crom_bind = lapply(
	crom, 
	function(x) rbind(
		x[[1]],
		x[[2]],
		x[[3]],
		x[[4]],
		x[[5]]))

crom_sample = lapply(
	crom_bind, function(x) split( x , f = x[[1]]))
	
crom_sample_IR = lapply(
	crom_sample, 
	function(x) lapply(
		x, 
		function(x) IRanges(x[[9]],x[[10]])))
	
crom_sample_Red = lapply(
	crom_sample_IR, 
	function(x) lapply(x, reduce))
	
crom_sample_dfH = lapply(
	crom_sample_Red, 
	function(x) lapply(x, data.frame))


#Para FROHs parciales

crom_sample_df12H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width < 2000000,]))
crom_sample_df24H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width >= 2000000 & x$width < 4000000, ]))
crom_sample_df48H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width >= 4000000 & x$width < 8000000, ]))
crom_sample_df816H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width >= 8000000 & x$width < 16000000, ]))
crom_sample_df16H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width >= 16000000, ]))


crom_sample_sum = lapply(crom_sample_dfH, function(x) lapply(x, function(x) sum(x[[3]])))
#En caso de querer estimar un FROH parcial XX hay que reemplazar "crom_sample_df" por "crom_sample_dfXX"
crom_sample_sum12 = lapply(crom_sample_df12H, function(x) lapply(x, function(x) sum(x[[3]])))
crom_sample_sum24 = lapply(crom_sample_df24H, function(x) lapply(x, function(x) sum(x[[3]])))
crom_sample_sum48 = lapply(crom_sample_df48H, function(x) lapply(x, function(x) sum(x[[3]])))
crom_sample_sum816 = lapply(crom_sample_df816H, function(x) lapply(x, function(x) sum(x[[3]])))
crom_sample_sum16 = lapply(crom_sample_df16H, function(x) lapply(x, function(x) sum(x[[3]])))

sumcrom = rbindlist(crom_sample_sum, fill=TRUE)
sumcrom[is.na(sumcrom)] <- 0
sumcrom = as.data.frame(sumcrom)
rownames(sumcrom) = c("01","10","11","12","13","14","15","16","17","18","19","02","20","21","22","23","24","25","26","27","28","29","03","04","05","06","07","08","09")
sumcrom = sumcrom[mixedsort(rownames(sumcrom)),mixedsort(colnames(sumcrom))]

sumcromH = sumcrom

#Repetir para los L
ff <- list.files(path="C:/Users/H/Desktop/RttPEDMAP/Seba_Espana/ROHs_N1_L",
                 full.names=TRUE, pattern = ".homozygousruns")
ff = ff[!ff %in% grep('500000', ff, value=T)] # Elimino los de 0.5 Mb
myfilelist <- lapply(ff, read.table, stringsAsFactors = F, sep = ",", header = T)
names(myfilelist) = unlist(lapply(ff, function(x) str_match(x,"[0-9]{1,2}_[0-9]{2}_[0-9]{7,8}")))
vector=rep(1:29, each=5)
crom = split(myfilelist,vector)
rm(vector)
crom_bind = lapply(crom, function(x) rbind(x[[1]],x[[2]],x[[3]],x[[4]],x[[5]]))
remove = c("RTT_11132","RTT_11133","RTT_11151","RTT_11161","RTT_11168","RTT_11172","RTT_11174") #Eliminadas por an?malas
crom_bind = lapply(crom_bind, function(x) x[! x$Label %in% remove, ])

#mismopa = c("RTT_11148","RTT_11134","RTT_11150","RTT_11129")
#crom_bind = lapply(crom_bind, function(x) x[x$Label %in% mismopa, ])

crom_sample = lapply(crom_bind, function(x) split( x , f = x[[1]]))
crom_sample_IR = lapply(crom_sample, function(x) lapply(x, function(x) IRanges(x[[9]],x[[10]])))
crom_sample_Red = lapply(crom_sample_IR, function(x) lapply(x, reduce))
crom_sample_dfL = lapply(crom_sample_Red, function(x) lapply(x, data.frame))


##################

#Para FROHs parciales

#crom_sample_df12 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width < 2000000,]))
#crom_sample_df24 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width >= 2000000 & x$width < 4000000, ]))
#crom_sample_df48 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width >= 4000000 & x$width < 8000000, ]))
#crom_sample_df816 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width >= 8000000 & x$width < 16000000, ]))
#crom_sample_df16 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width >= 16000000, ]))


crom_sample_sum = lapply(crom_sample_dfL, function(x) lapply(x, function(x) sum(x[[3]])))
#En caso de querer estimar un FROH parcial XX hay que reemplazar "crom_sample_df" por "crom_sample_dfXX"
sumcrom = rbindlist(crom_sample_sum, fill=TRUE)
sumcrom[is.na(sumcrom)] <- 0
sumcrom = as.data.frame(sumcrom)
rownames(sumcrom) = c("01","10","11","12","13","14","15","16","17","18","19","02","20","21","22","23","24","25","26","27","28","29","03","04","05","06","07","08","09")
sumcrom = sumcrom[mixedsort(rownames(sumcrom)),mixedsort(colnames(sumcrom))]

sumcromL = sumcrom

sumcrom = cbind(sumcromH, sumcromL)#H y L grupos de baja y alta endogamia(no usar en cabras)

largos = seqlengths(BSgenome.Btaurus.UCSC.bosTau8)
FROHSxChr = sapply(sumcrom, function(x) x/largos[1:29])
FROHSxChr = t(FROHSxChr)
FROHSxChr = as.data.frame(FROHSxChr)
colnames(FROHSxChr) = c(1:29)
FROHSxChr$Group = c(rep("HI",ncol(sumcromH)),rep("LI",ncol(sumcromL)))
FROHSxChr$Sample = rownames(FROHSxChr)

FROHSxChr = melt(FROHSxChr, id.vars=c("Group", "Sample"))
colnames(FROHSxChr) = c("Group","Sample","Chromosome","FROH")


#Graficar FROH por cromosoma
ggplot(FROHSxChr, aes(x=Chromosome, y=FROH, fill=Group)) + theme_classic() + 
  geom_boxplot() + ylab(expression(paste("F"[ROH]))) + 
  stat_summary(fun.y="mean", geom="point", size=1.5,
               position=position_dodge(width=0.75), color="white") + 
  scale_x_discrete(name ="Chromosome")

#Para obtener los FROH de la longitud que hayamos especificado en crom_sample_df"XX"
sums = colSums(sumcrom)
largo = sum(as.numeric(largos [1:29]))
FROHs = sums/largo
FROHs = FROHs[order(names(FROHs))]
FROHs