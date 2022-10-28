library(tidyverse)
library(ggfortify)



##CHANGE ANNOTATION TO BASELINE MAKE IT CONSISTENT VISIT 1, VISIT 2

a<-read_tsv("All_Proteins_groups2B.txt")
#a<-read_tsv("All_Proteins_groups2Minus_10.1.txt")

a<-a %>% unite("z", Sample.ID:TimePoint, remove = FALSE)

a$TimePoint<-as.factor(a$TimePoint)

#a<-a %>%  mutate(TimePoint, TimePoint = ifelse(Group == 'AMC', "AMC", TimePoint))

a<-a %>%  mutate(TimePoint2= case_when(
   (Group == 'AMC' ~ "Healthy Controls"),
   (TimePoint == '0' ~ "Patients (Baseline)"),
   
   (TimePoint == '1' ~ "Patients (Days 8-15)"),
   (TimePoint == '2' ~ "Patients (Days 16-30)"),
   (TimePoint == '3' ~ "Patients (Days 48-59)")))
  
  
  

#a<-a %>% filter(Group=="PT")%>% filter(TimePoint != "NA")
a<-a %>% filter(Group %in% c("PT","AMC"))%>% filter(TimePoint2 != "NA")


###REMOVE 10.1
#a<-a %>% filter(!str_detect(Sample.ID, 'P10.1 PT'))

rownames(a)<-make.names(a$z,unique = TRUE)

#prcomp(a[7:length(a)])

pc<-prcomp(log2(a[8:length(a)-1]))
#pca$x[,1:2]
res.pca<-pc

geom_line(data=pc,aes(x=pc$x[,1],y=pc$x[,2]))

pall<-autoplot(pc,scale=TRUE,text=a[1] , shape = TRUE,
               main="Total 94 proteins different at baseline", #PT + AMC",
         data=a,
         label = FALSE,
         label.size = 3,
         text = paste("ID:", a$Sample.ID),
         colour="TimePoint2",frame=T,label.show.legend = FALSE
) 




#################################################


grp<-read_tsv("group_test")


down<-grp %>% filter(`Fold Change` < 0 ) %>% select(id)

up<-grp %>% filter(`Fold Change` > 0 ) %>% select(id)

downdat<-a %>% select(1:6,TimePoint2,dput(down$id))

pc<-prcomp(log2(downdat[8:length(downdat)]))
#pca$x[,1:2]
res.pca<-pc

geom_line(data=pc,aes(x=pc$x[,1],y=pc$x[,2]))

pl<-autoplot(pc,scale=TRUE,
             #text=downdat[1] ,
             shape = TRUE,
         data=downdat,main="26 proteins down at baseline", ##PT Vs AMC = Minus FC",
         label = FALSE,
         #label.size = 3,
        # text = paste("ID:", downdat$Sample.ID),
         colour="TimePoint2",frame=T,legend = FALSE,label.show.legend = FALSE
)
updat<-a %>% select(1:6,TimePoint2,dput(up$id))



pc<-prcomp(log2(updat[8:length(updat)]))
#pca$x[,1:2]
res.pca<-pc

geom_line(data=pc,aes(x=pc$x[,1],y=pc$x[,2]))

ph<-autoplot(pc,scale=TRUE,text=updat[1] , shape = TRUE,
         data=updat,main="68 proteins up at baseline",## PT Vs AMC = Plus FC",
         label = FALSE,
         label.size = 3,
         text = paste("ID:", updat$Sample.ID),
         colour="TimePoint2",frame=T,legend = FALSE,label.show.legend = FALSE
)


#pdf("pca_minus_10_1.pdf", width = 11)
pdf("pca2.pdf", width = 11)
library("gridExtra")
#grid.arrange(pall + theme(legend.title = element_blank()) +   scale_color_discrete(breaks=c("AMC","Patients (Baseline)","Patients (Days 8-15)", "Patients (Days 48-59)")),                             # First row with one plot spaning over 2 columns
 ##            arrangeGrob(pl+ theme(legend.title = element_blank()) +   scale_color_discrete(breaks=c("AMC","Patients (Baseline)","Patients (Days 8-15)", "Patients (Days 48-59)")) ,
#                         ph+ theme(legend.title = element_blank()) +   scale_color_discrete(breaks=c("AMC","Patients (Baseline)","Patients (Days 8-15)", "Patients (Days 48-59)")), ncol = 2), # Second row with 2 plots in 2 different columns
#             nrow = 2) 


grid.arrange(pall + theme(legend.title  = element_blank()),                             # First row with one plot spaning over 2 columns
            arrangeGrob(pl+   theme(legend.title  = element_blank())   ,ph+   theme(legend.title  = element_blank())  , 
                        ncol = 2), # Second row with 2 plots in 2 different columns
             nrow = 2) 



dev.off()





###grid.arrange(pall + theme(legend.title = element_blank()) +   scale_color_discrete(breaks=c("AMC","Patients (Baseline)","Patients (Days 8-15)", "Patients (Days 48-59)")),                             # First row with one plot spaning over 2 columns
##            arrangeGrob(pl+ theme(legend.title = element_blank()) +   scale_color_discrete(breaks=c("AMC","Patients (Baseline)","Patients (Days 8-15)", "Patients (Days 48-59)")) ,
#                         ph+ theme(legend.title = element_blank()) +   scale_color_discrete(breaks=c("AMC","Patients (Baseline)","Patients (Days 8-15)", "Patients (Days 48-59)")), ncol = 2), # Second row with 2 plots in 2 different columns
#             nrow = 2) 



####################################

###WRONG IDS
pc<-prcomp(log2(a[8:length(a)-1]))
df_out <- as.data.frame(pc$x)


df_out$Time<-sapply( strsplit(as.character(row.names(a)), "_"), "[[", 2 )
df_out$ID<-sapply( strsplit(as.character(row.names(a)), "_"), "[[", 1 )
##df_out$ID2<-as.character(row.names(a))
ggplot(df_out,aes(x=PC1,y=PC2,color=Time,label=ID,group=ID ))+geom_point() + 
     geom_text() +
     geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) 

