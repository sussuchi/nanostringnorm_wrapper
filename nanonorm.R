######### handling configs, libraries and inputs ###########

## to run wrapped.R: Rscript --default-packages=methods,utils ${__tool_directory__}/nanonorm.R --rccfile "$i_rccfile" --output $o_output --outpdf $o_outpdf --codecount $codecount --background $background --samplecont $samplecont --othernorm $othernorm --boxplot $boxplot --controls $controls

## Set up R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

## Avoid crashing Galaxy with an UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Import required libraries
#suppressPackageStartupMessages({
library(WriteXLS)
library(getopt)
library(gdata)
library(psych)
library(ComplexHeatmap)
library(circlize)
library(XML)
library(RColorBrewer)
library(NanoStringNorm)
#])
## Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

## Get options using the spec as defined by the enclosed list
## Read the options from the default: commandArgs(TRUE)
option_specification = matrix(c(
  'rccfile', 'i', 1, 'character',
  'output', 'n', 2, 'character',
  'outpdf', 'p', 2, 'character',
  'codecount', 'c', 1, 'character',
  'background', 'b', 1, 'character',
  'samplecont', 's', 1, 'character',
  'othernorm', 'o', 1, 'character',
  'boxplot', 'l', 1, 'character',
  'controls', 't', 1, 'character',
  'househk', 'h', 1, 'character',
  'trait', 'tr', 1, 'character',
  'statistic', 'st', 1, 'character'
  
), byrow=TRUE, ncol=4)


## Parse options
options <- getopt(option_specification)

## Ler arquivo xlsx
rawdata <- read.xls.RCC(options$rccfile, sheet =1)
data <- rawdata$x

## ler tabela de traits e reordenar
class <- read.xls(options$trait, sheet = 1)
rownames(class) <- class[,1]
traits <- class[,-1]
samplenames <- names(data)[-c(1:3)]
traits.final <- traits[match(samplenames,rownames(traits)),]

## plot dos controles
datacontrols <- subset(data, Code.Class %in% c("Positive","Negative"))
datacontrols$color[datacontrols$Code.Class=="Positive"] <- "red"
datacontrols$color[datacontrols$Code.Class=="Negative"] <- "blue"
n <- ncol(datacontrols)-1
pdf(options$controls, width = 10, height = 90)
dotchart(log2(as.matrix(datacontrols[4:n])),labels = datacontrols$Name, cex = .7, main ="Controls before normalization", xlab = "log2 of raw counts", gcolor = "black", color = datacontrols$color)
dev.off()


# normalização
data.norm.fs <- NanoStringNorm(
  x = data,
  CodeCount = options$codecount,
  Background = options$background,
  SampleContent = options$samplecont,
  OtherNorm = options$othernorm,
  round.values = TRUE,
  take.log = TRUE,
  traits = traits.final
)

WriteXLS(data.norm.fs$normalized.data, options$output)
statistic <- data.norm.fs$gene.summary.stats.norm
# linha de fix da tabela
statistic <- data.frame(Genes = row.names(statistic), statistic)
WriteXLS(statistic, options$statistic)



# boxplot
matrix <- data.norm.fs$normalized.data
dataEndogenous <- subset(matrix, Code.Class %in% "Endogenous")
dataEndogenous <- dataEndogenous[,-1:-3]
matrixbefore <- data
rownames(matrix) <- data$Name
dataEndogenousbefore <- subset(matrixbefore, Code.Class %in% "Endogenous")
dataEndogenousbefore <- dataEndogenousbefore[,-1:-3]
dataEndogenousbefore <- log2(dataEndogenousbefore)
pdf(options$boxplot, width = 14);
boxplot(dataEndogenous, col = "blue", las=2, main = "Endogenous after normalization");
boxplot(dataEndogenousbefore, col = "blue", las=2, main = "Endogenous before normalization");
dev.off()


# boxplot hk
dataHK <- subset(matrix, Code.Class %in% "Housekeeping")
dataHK <- dataHK[,-1:-3]
dataHK <- t(dataHK)
dataHKbefore <- data
rownames(dataHKbefore) <- data$Name
dataHKbefore <- subset(dataHKbefore, Code.Class %in% "Housekeeping")
dataHKbefore <- dataHKbefore[,-1:-3]
dataHKbefore <- t(log2(dataHKbefore))
pdf(options$househk, width = 14);
boxplot(dataHK, col = "red", las=2, main = "Housekeeping after normalization");
boxplot(dataHKbefore, col = "red", las=2, main = "Housekeeping before normalization");
dev.off()

# plot all nanostringnorm 
pdf(options$outpdf);
Plot.NanoStringNorm(
  x = data.norm.fs,
  label.best.guess = TRUE,
  plot.type = 'all'
);
dev.off()

