
###CHECK THE NUMBERS OF BASELINE VERSUS V1 LOOKS TOO MANY
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
a<-read_tsv("chapel_correct.txt")

           


library(readxl)
data <- read_excel("paired_max_visit_PT.xlsx")
#data <- read_excel("A Generic Plotting and Statistics Web Tool Version 0.01-4.xlsx")
#data <- read_excel("paired_baseline_V1.xlsx")

data<-data %>% filter(`Wilcoxin Rank Sum Test (fdr)` < 0.05)


a$TimePoint<-as.numeric(as.character(a$TimePoint))
b<-a %>% mutate(visit= case_when(
  # (TimePoint <= 0  ~ "0"),(TimePoint == 9 ~ "1"),(TimePoint == 8 ~ "1"),
  # (TimePoint <= 10 & TimePoint > 0 ~ "1"),
  # (TimePoint <= 20 & TimePoint > 10 ~ "2"),
  # (TimePoint <= 30 & TimePoint > 20 ~ "3"),
  # (TimePoint <= 40 & TimePoint > 30 ~ "4"),
  # (TimePoint <= 50  & TimePoint > 40 ~ "5"),
  # (TimePoint <= 60  & TimePoint > 50 ~ "6"), (TRUE  ~ "0")
  
  (TimePoint <= 0  ~ "0"),
  (TimePoint <= 15 & TimePoint >= 8 ~ "1"),
  (TimePoint <= 30  & TimePoint >= 16 ~ "2"),
  (TimePoint <= 59  & TimePoint >= 48 ~ "3"),
  str_detect(SampleId,'AMC') ~ "0",
  (TRUE  ~ "NA")
  
  
  
) ) %>%
  mutate(Subjects = str_replace(SampleId, " D.*", "")) %>%
  mutate(SampleGroup= str_replace_all(SampleId, regex(".*(AMC|PT|PC).*"), "\\1"  ) )
#c<-b %>% filter(SampleNotes != 'No sample' | is.na(SampleNotes) )%>% dplyr::select(Subjects,visit,SampleGroup,Age,Gender,
#                                                                                   "EGFR","ENG","NTRK2","CNTN1","APOL1","NOTCH1","CNTFR","IL1R2",
#                                                                                   "ALCAM","BMP1","RET","CD86","CD109","L1CAM","NTRK3","ADGRE2","CHL1",
#                                                                                   "REG4","SELL","MMP3","SET","GSK3A.GSK3B","A2M","FGG","TIMP1","APOM",
#                                                                                   "IL2","FGA.FGB.FGG.1","IGHG1","NXPH1","GSN","PAPPA","CA6","IGFBP2",
#                                                                                   "FABP3")

c<-b %>% filter(SampleNotes != 'No sample' | is.na(SampleNotes) ) %>% dplyr::select(Subjects,visit,SampleGroup,Age,Gender,SampleNotes, data$id)

c<- c %>% dplyr::filter(!(Subjects == 'P1.1 PC' & visit == 0) )

#Code to take top/oldest visit parsed from a rqnge of visits
#d<-c  %>%group_by(Subjects,visit) %>%slice(which.max(visit))  %>% arrange(desc(visit)) %>%  group_by(Subjects) %>% top_n(1,visit) %>% filter(visit > 0)


ddd<-c %>% filter(visit ==3)
d<-c %>% filter(visit ==1)
dd<-c %>% filter(visit ==2)
e<-c %>% filter(visit ==0)


f<-merge(d,e, by=colnames(d), all=TRUE)[colnames(d)]
ff<-merge(dd,f, by=colnames(dd), all=TRUE)[colnames(dd)]
#ADD #RD TIMEPOINT HERE THEN MERGED
fff<-merge(ddd,ff, by=colnames(ddd), all=TRUE)[colnames(ddd)]

f<-fff
##c<- c%>% dplyr::filter(SampleGroup %in% c("AMC","PT"))
#f<- f%>% dplyr::filter(SampleGroup %in% c("AMC","PT","PC"))
f<- f%>% dplyr::filter(SampleGroup %in% c("AMC","PT"))
#c<-c %>% dplyr::filter(visit %in% c(0,2))
#f<-f %>% dplyr::filter(visit %in% c(0,2))
f<- f %>% dplyr::arrange(SampleGroup,visit)
#b<- b%>% filter(SampleGroup %in% c("AMC","PT"))
#b <- b %>% filter(TimePoint == 0)
ha_rot_cn = HeatmapAnnotation(text = anno_text(f$Subjects, rot = 90))
#offset = unit(-0.5, "cm")))

#ha_rot_cn = HeatmapAnnotation(text = anno_text(c$Subjects, rot = 90))

#ha <- HeatmapAnnotation(df = data.frame(Group = c$SampleGroup),visit=c$visit)
ha <- HeatmapAnnotation(df = data.frame(Group = f$SampleGroup), Gender=f$Gender,Visit=f$visit, 
                        #SampleNotes=f$SampleNotes,
                        Age=anno_barplot(f$Age),col=col )

boxplotCol <- HeatmapAnnotation(boxplot=anno_boxplot(data.matrix(t(1-scale(f[7:length(data$id)])* -1)), border=TRUE, gp=gpar(fill="#CCCCCC"), pch=".", size=unit(2, "mm"), axis=TRUE, annotation_width=unit(c(1, 5.0), "cm"), which="col"))

##hea1<-Heatmap(t(1-scale(f[6:length(data$id)])* -1), top_annotation = ha,bottom_annotation= ha_rot_cn,show_row_names = TRUE, column_title = 'PT: Visit 0 Vs Last Visit ',
##           cluster_columns =  cluster_within_group(t(f[6:length(data$id)]),f$SampleGroup))
f$Subjects<-factor(f$Subjects,levels=c("P1.1 AMC", "P10.1 AMC", "P13.1 AMC", "P13.2 AMC", "P16.1 AMC", 
                                      "P17.1 AMC", "P2.1 AMC", "P9.1 AMC", "P1.1 PT", "P10.1 PT", "P13.1 PT", 
                                       "P13.2 PT", "P16.1 PT", "P17.1 PT", "P2.1 PT", "P9.1 PT"))


f<-f %>% rename(IGHG1 = IGHG1.IGHG2.IGHG3.IGHG4.IGK..IGL.)

hea1<-Heatmap(t(1-scale(f[7:length(f)])* -1), heatmap_legend_param = list(
  title = "Z-score",at = c(4,3,2,1,0,-1),legend_height = unit(3, "cm"),color_bar = "continuous" ),
     cluster_column_slices =FALSE,
       cluster_columns = FALSE,#col=col_fun,
  column_split = rep(c( "Healthy Controls","Patients (Baseline)","Patients (Days 08-15)", "Patients (Days 16-30)","Patients (Days 48-59)"), c(8,8,8,7,8)),
  # column_split = rep(c( "Healthy Controls","Patients (Baseline)", "Patients (Visit2)","Patients (Visit3)","Patients (Visit4)"), c(8,8,8,7,8)),
        row_split=2,show_heatmap_legend = TRUE,name="Value",border = TRUE, show_column_dend = FALSE,
        show_row_dend = FALSE, row_names_side  = "left",column_title_side = "bottom",row_title = NULL,
        height = unit(0.5, "cm")*ncol(f),width = unit(0.7, "cm")*nrow(f))
#, row_title = NULL
# column_labels = rep(c( "Healthy Controls","Baseline", "Days 8-15","Days 48-59"), c(8,8,8,8)),


#hea1<-Heatmap(t(1-scale(f[7:length(f)])* -1), 
#              +       #  bottom_annotation= ha_rot_cn,
#                +      #   show_row_names = TRUE,
#                +        # column_title = ' ',
#                +      # cluster_columns = FALSE,column_split = rep(c( "A","B","C","D"), c(8,8,8,8)),
#                +      show_column_names = TRUE,cluster_column_slices =FALSE,
##              +        cluster_columns = FALSE,column_dend_reorder = FALSE,
#              +    #  column_split = rep(c( " "," ", " "," "), c(8,8,8,8)),
#                +      column_split = rep(c( "Healthy Controls","Baseline", "Days 8-15","Days 48-59"), c(8,8,8,8)),
#              +     # column_labels = rep(c( "Healthy Controls","Baseline", "Days 8-15","Days 48-59"), c(8,8,8,8)),
#                +         row_split=2,show_heatmap_legend = TRUE,name="Value",border = TRUE, show_column_dend = FALSE,
#              +         show_row_dend = FALSE, row_names_side  = "left",column_title_side = "bottom",
#              +         height = unit(0.5, "cm")*ncol(f),width = unit(0.7, "cm")*nrow(f))
#>



####################################################################################


library(readxl)
data <- read_tsv("group_test")
#data <- read_excel("D0_AMC_PT_0.1_SD.xlsx")



library(ComplexHeatmap)
library(circlize)
library(tidyverse)
a<-read_tsv("chapel_correct.txt")




a$TimePoint<-as.numeric(as.character(a$TimePoint))
b<-a %>% mutate(visit= case_when(
  (TimePoint <= 0  ~ "0"),(TimePoint == 9 ~ "1"),(TimePoint == 8 ~ "1"),
  (TimePoint <= 10 & TimePoint > 0 ~ "1"),
  (TimePoint <= 20 & TimePoint > 10 ~ "2"),
  (TimePoint <= 30 & TimePoint > 20 ~ "3"),
  (TimePoint <= 40 & TimePoint > 30 ~ "4"),
  (TimePoint <= 50  & TimePoint > 40 ~ "5"),
  (TimePoint <= 60  & TimePoint > 50 ~ "6"), (TRUE  ~ "0")
) ) %>%
  mutate(Subjects = str_replace(SampleId, " D.*", "")) %>%
  mutate(SampleGroup= str_replace_all(SampleId, regex(".*(AMC|PT|PC).*"), "\\1"  ) )
##c<-b %>% dplyr::select(Subjects,visit,SampleGroup,REG4,APOL1,APOM,IL1R2,ANXA5,TIMP1,KIT,NTRK3,EGFR,PROS1,CNTFR,C1QA,CYP3A4,PROC,CNTN1,RET,CNTN4,BMP1,CA6,CTSV,SELE,ECE1,A2M,PSMA1,ING1,KLKB1,CFH,L1CAM,ADGRE2,PLXNC1,GSN,SELL,NTRK2,ENG,NOTCH1,NRCAM,TYK2,CD86,SET,TNNI3,ALCAM,EFEMP1,CRYZL1,IL6ST,IGHG1,LSAMP,CCDC80,KDR,IGF2R,SKP1,COL8A1,PARK7,CCL11,NRP1,IGFBP2,NEGR1,GHR,WFIKKN1,CD109,CST6,HIPK3,CHL1,NACA,NCAM1,ERP29,GPNMB,CD200R1,SLAMF6,FGFR1,EHMT2,LRP1B)
###c<-b %>% dplyr::select(Subjects,visit,SampleGroup,"EGFR","ENG","NTRK2","CNTN1","APOL1","NOTCH1","CNTFR","IL1R2","ALCAM","BMP1","RET","CD86","CD109","L1CAM","NTRK3","ADGRE2","CHL1","REG4","SELL","MMP3","SET","GSK3A.GSK3B","A2M","FGG","TIMP1","APOM","IL2","FGA.FGB.FGG.1","IGHG1","NXPH1","GSN","PAPPA","CA6","IGFBP2","FABP3")

c<-b %>% filter(SampleNotes != 'No sample' | is.na(SampleNotes) ) %>% dplyr::select(Subjects,visit,SampleGroup,Gender,Age,SampleNotes,data$id)

c<- c %>% dplyr::filter(!(Subjects == 'P1.1 PC' & visit == 0) )

c<- c%>% dplyr::filter(SampleGroup %in% c("PT","AMC")) #, "PC"))

#c<-c %>% dplyr::filter(visit %in% c(0,2))
c<-c %>% dplyr::filter(visit %in% c(0))
c<- c %>% dplyr::arrange(SampleGroup,visit)
#b<- b%>% filter(SampleGroup %in% c("AMC","PT"))
#b <- b %>% filter(TimePoint == 0)
c<-c %>% rename(IGHG1 = IGHG1.IGHG2.IGHG3.IGHG4.IGK..IGL.)
#c<-c %>% rename(IGH1 = IGHG1.IGHG2.IGHG3.IGHG4.IGK..IGL.)

ha_rot_cn = HeatmapAnnotation(text = anno_text(c$Subjects, rot = 90))

ha <- HeatmapAnnotation(df = data.frame(Group = c$SampleGroup),Gender=c$Gender,Visit=c$visit ,
                        #SampleNotes=c$SampleNotes,
                        col=col, Age=anno_barplot(c$Age),show_legend = FALSE)
#hea<-Heatmap(t(1-scale(c[4:38])* -1), top_annotation = ha,bottom_annotation= ha_rot_cn,show_row_names = TRUE, column_title = 'PT vs PC: Baseline', cluster_columns = TRUE

#f$Subjects<-factor(f$Subjects,levels=c("P1.1 AMC", "P10.1 AMC", "P13.1 AMC", "P13.2 AMC", "P16.1 AMC", 
#"P17.1 AMC", "P2.1 AMC", "P9.1 AMC", "P1.1 PT", "P10.1 PT", "P13.1 PT", 
#"P13.2 PT", "P16.1 PT", "P17.1 PT", "P2.1 PT", "P9.1 PT"))



hea<-Heatmap(t(1-scale(c[7:length(c)])* -1), heatmap_legend_param = list(
  title = "Z-score",at = c(4,3,2,1,0,-1),legend_height = unit(3, "cm"),color_bar = "continuous" ),
          #   top_annotation = ha,
          #   bottom_annotation= ha_rot_cn,
             show_row_names = TRUE,# column_title = 'Baseline',
             cluster_columns = TRUE,#col=col_fun, 
          #   column_split = rep(c("C", "D", "E"), c(8,7,8)), row_split=2, show_heatmap_legend = TRUE,name="Value",border = TRUE,
          column_split = rep(c("Healthy Controls", "Patients"), c(8,8)), row_split=2, show_heatmap_legend = TRUE,name="Value",border = TRUE,
              show_column_dend =FALSE,column_title_side = "bottom",
             show_row_dend = FALSE,
             height = unit(0.5, "cm")*ncol(c),width = unit(0.7, "cm")*nrow(c),
            row_names_side  = "left",row_names_gp =gpar(just="left"),row_title = NULL
          #,rowAnnotation(just = "left")
            )
             #heatmap_height = unit(0.5, "cm")*ncol(c),heatmap_width = unit(0.7, "cm")*nrow(c))
             
####ha %v% NULL

#pdf("all.pdf",height=15)

#draw(hea1 + hea, padding = unit(c(20, 20, 10, 10), "mm"),column_title = "?")
#dev.off()




library(gridExtra)
pdf("test2_v2.pdf",height=23, width=27)
#png("test2_v2.png",height=2300, width=2700)
#tiff("test2_v2.tiff",height=2300, width=2700)
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2,widths  = unit(c(3.5,7.5), "null"),just = "top")))
#pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2,heights = unit(c(0.5,10), "null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1,just = "top"))

draw(hea, newpage = FALSE)
#decorate_heatmap_body("aaa", {
  
 # grid.rect (c(0.5,1.52), c(-0.73,-0.73),height=3.45,gp = gpar(fill = "transparent", col = "black", lty=2,lwd = 2)  )
  
#})
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,just = "top"))

##draw(hea1, newpage = FALSE)
draw(hea1, newpage = FALSE)

#decorate_heatmap_body("hhh", {
  
 # grid.rect (c(0.5,1.52), c(0.11,0.11),height=1.78,gp = gpar(fill = "transparent", col = "black", lty=2,lwd = 2)  )
  
#})
upViewport()
#grid.text("Chaple Protocol: SomaLogic Results", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2) )
dev.off()

##f$Subjects<-factor(f$Subjects,levels=c("P1.1 AMC", "P10.1 AMC", "P13.1 AMC", "P13.2 AMC", "P16.1 AMC", 
"P17.1 AMC", "P2.1 AMC", "P9.1 AMC", "P1.1 PT", "P10.1 PT", "P13.1 PT", 
"P13.2 PT", "P16.1 PT", "P17.1 PT", "P2.1 PT", "P9.1 PT"))



################################################

Heatmap(t(1-scale(f[7:length(data$id)])* -1), top_annotation = ha,bottom_annotation= ha_rot_cn,show_row_names = TRUE, column_title = 'Baseline Vs Last Visit',
        cluster_columns = TRUE,column_split = rep(c("C", "D", "E","F", "G"), c(8,7,8,8,8)), row_split=2,show_heatmap_legend = FALSE,name="hhh")

decorate_heatmap_body("hhh", {
  
  grid.rect (c(0.5,4.85), c(0.1,0.1),height=1.8,gp = gpar(fill = "transparent", col = "black", lwd = 1.5)  )
  
})



he<-  Heatmap(t(1-scale(c[7:length(data$id)])* -1),# top_annotation = ha,bottom_annotation= ha_rot_cn,
              show_row_names = TRUE,column_title = 'Baseline Differences PT Vs AMC', cluster_columns = TRUE,
        column_split = rep(c("C", "D", "E"), c(8,7,8)), row_split=2, show_heatmap_legend = FALSE,name="aaa") 
  he 
  decorate_heatmap_body("aaa", {
          
                grid.rect (c(0.5,1.55), c(-0.75,-0.75),height=3.5,gp = gpar(fill = "transparent", col = "black", lwd = 1.5)  )
           
            })



  
  
  
  ###################
  
 ### HeatmapAnnotation(df = data.frame(Group = f$SampleGroup), Gender=f$Gender,Visit=f$visit,Age=anno_barplot(f$Age), SampleNotes=f$SampleNotes,text = anno_text(f$Subjects, rot = 90)) %v% NULL