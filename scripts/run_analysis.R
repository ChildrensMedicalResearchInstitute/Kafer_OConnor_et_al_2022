#! /usr/bin/env Rscript

# Author(s): Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


#Test for standard library dependencies

list.of.packages <- c("ggplot2","tidyverse","reshape2","combinat","multcomp","MASS","boot","forestplot","haven")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cran.csiro.au/')

#load packages
library(ggplot2)
library(tidyverse)
library(reshape2)
library(multcomp)
library(MASS)
library(boot)
library(forestplot)
library(zoo)

# Upload files to the data directory in comma separated format (excel save as...)
base_dir <- here::here()
data_file <- file.path( base_dir, "data", "Chromosome_Spreads.csv.zip")
print(data_file)
df.raw<-readr::read_csv(file=data_file, col_names=TRUE) #load data to a data frame (df)
df.raw$Sample <- factor(df.raw$Sample, levels = c("Control","0.25 mM_HU","2 mM_HU"))
df.raw$Cell_Line <- factor(df.raw$Cell_Line, levels = c("ES", "Epi"))
df.raw$Phenotype <- factor(df.raw$Phenotype, levels = c("Breaks", "Separated_Fragmented", "Normal mitosis"))

df.totals <- df.raw %>%
group_by(Repeat,Cell_Line,Sample) %>%
summarise(total=sum(Count))

#%%
df <- merge(df.raw, df.totals,by=c("Repeat","Cell_Line","Sample")) %>%
mutate(Norm=Count/total)

# multivariate summary
df.summary2 <- df %>%
  group_by(Cell_Line,Phenotype,Sample) %>%  # group by two variables: Sample and Class
  summarise(
    n=n(),
    mean=mean(Norm),
    sd=sd(Norm)
  )  %>%
  mutate( se=sd/sqrt(n))

df.summary2$mean2=df.summary2$mean # add new column mean2
df.summary2[df.summary2$Cell_Line=="ES","mean2"]=-df.summary2[df.summary2$Cell_Line=="ES","mean2"] #change sign of ES

figure_file <- file.path( base_dir, "figures", "stackedbarchart.pdf")
print(figure_file)
pdf(figure_file, width=12, height=5)

options(repr.plot.width=7, repr.plot.height=3)

# Here we use standard error for uncertanty
ggplot(df.summary2,aes(x=Sample, y=mean, fill=Phenotype)) + # fill defines the variable for subgroups
  geom_bar(stat="identity", color="black") +
  ggtitle("") +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ #removes background
scale_fill_manual(values=c("#800000", "#006400", "#008B8B","#4682B4","#D2691E"))+
facet_wrap(~Cell_Line, scales="fixed")+
labs(fill = "Phenotype", y="Mean proportion",x="Sample")

df.cast <- dcast(df, Repeat+Cell_Line + Phenotype ~ Sample, value.var="Norm")

calculate_difference <- function(ds,rows,columns,denom){
    #This function calculate the ratio between each column and a given denominator (column)
    ds.numerator<-ds[,columns]
    ds.denom<-ds[,denom]
    nc<-ncol(ds.numerator) #number of columns
    nr<-nrow(ds) #number of rows
    result<-data.frame(Sets=seq(1,nr**2)) #the result has nr**2 (combination of every column divided by every value of the denominator)
  result<-cbind(result,rows)
  for(i in 1:nc)  # iterate every column except the first (expected index in first column)
    {
        data<-c() # empty vector to store the result
        for(j in 1:nr){ #iterate each row of the numerator
        for(k in 1:nr){ #iterate each row of the denominator
            data<-rbind(data,ds.numerator[j,i]-ds.denom[k,]) # calculate the log 2 ratio
        }
        }
        result<-cbind(result, data) # add a new column with the values
    }
    return (result) #return the result
}
df.diff<-df.cast %>%
  group_by(Phenotype,Cell_Line) %>%  # group by two variables: Pahse and Cell_Line
  group_map(calculate_difference,columns=c("0.25 mM_HU","2 mM_HU"),denom="Control") %>%  # Apply calculate_diff to each group with "Control as reference"
  melt(id.var = c('Cell_Line','Phenotype','Sets'), variable.name = 'Sample') # reshape the result to have columns Cell_Line, Phase, Sets, Sample and value


round_p_value<-function(p_values,min_value){
    result<-c()
    for( p in p_values){
        if(is.na(p)){
            result<-rbind(result,"")
            next
        }
        else{
        if(p < min_value){
            result<-rbind(result,paste("p < ",min_value))
        }
        else{
            result<-rbind(result,as.character(round(p,digits = abs(round(log10(min_value))))))
        }
        }
    }
    return(result)
}


# multivariate summary
df.summary3 <- df.diff %>%
  group_by(Cell_Line,Phenotype,Sample) %>%  # group by two variables: Sample and Class
  summarise(
    n=n(),
    mean=mean(value),
    sd=sd(value),
    ci_l=abc.ci(value, weighted.mean)[2],
    ci_u=abc.ci(value, weighted.mean)[3],
    p=wilcox.test(value,exact=FALSE)$p.value
  )  %>%
  mutate( se=sd/sqrt(n))

df.summary3[,"p_adjust"]=p.adjust(df.summary3$p,method = "BH")  # c("holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini & Hochberg), "BY","fdr", "none")
df.summary3[,"p_adjust_fmt"]=round_p_value(df.summary3$p_adjust,min_value=0.001)
df.summary3$ci_l<-na.fill(df.summary3$ci_l,fill=0)
df.summary3$ci_u<-na.fill(df.summary3$ci_u,fill=0)


figure_file <- file.path( base_dir, "figures", "forestplot.pdf")
print(figure_file)

pdf(figure_file, width=5, height=3)

format_column<-function(values){
    len_values<-length(values)
    result<-as.vector(rep(NA,len_values))
    current<-""
    for(i in 1:len_values){
       value_str<-toString(values[i])
        if(current != value_str){
            current = value_str
            result[i]=current
        }
    }
    return (result)
}

options(repr.plot.width=8, repr.plot.height=3)
#make a data frame with mean, lower and upper 95% confidence interval bounds
forestplot_data<-cbind(
    c(NA,as.vector(df.summary3$mean)),
    c(NA,as.vector(df.summary3$ci_l)),
    c(NA,as.vector(df.summary3$ci_u))
     )

#make table with additional information
tabletext<-cbind(
    c("Cell Line",format_column(df.summary3$Cell_Line))
    ,c("Phenotype",format_column(df.summary3$Phenotype))
    ,c("Sample",as.vector(df.summary3$Sample))
    ,c("Mean",as.vector(round(df.summary3$mean,digits = 2)))
    ,c("CI",
       paste("["
        ,c(as.vector(round(df.summary3$ci_l,digits = 2)))
             ,", "
        ,c(as.vector(round(df.summary3$ci_u,digits = 2)))
        ,"]")
       )
    ,c("p-value",as.vector(df.summary3$p_adjust_fmt))
    ,c("n",as.vector(df.summary3$n))
     )
#create a forest plot
forestplot(
    tabletext
    ,forestplot_data
    ,new_page = FALSE
    ,hrzl_lines = list("2"=gpar(col="#444444",lwd=2)
                             ,"4"=gpar(col="#444444",columns=2:8)
                             ,"6"=gpar(col="#444444",columns=2:8)
                             ,"8"=gpar(col="#444444",lwd=2)
                             ,"10"=gpar(col="#444444",columns=2:8)
                             ,"12"=gpar(col="#444444",columns=2:8)
                             ,"14"=gpar(col="#444444",lwd=2)
                            ),
           txt_gp = fpTxtGp(label = gpar(fontsize=6),summary=gpar(fontsize=8),xlab=gpar(fontsize=12)),
           graph.pos = 4,
           line.margin = .1, # We need to add this to avoid crowding
           lwd.zero=2,
           #is.summary=c(as.vector(df_s4$summary)),
            boxsize = .125,
           clip=c(-Inf,Inf),
           colgap=unit(2,"mm"),
           xlog=FALSE,
           grid=TRUE,
           xlab="Mean difference",
            vertices = TRUE,
           graphwidth=unit(15,"mm"),
           col=fpColors(box="#4682B4",line="#4682B4", summary="#B22222")
          )



