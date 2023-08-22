#loading libraries to work with

library(readr)
library(readxl)  # to import Microsoft Excel (.xlsx) files
library(dplyr)
library(stringr)
library(tidyverse) # data manipulation and plotting packages
library(scales)
library(gtools) # used to sort `well` column
library(viridis)
library(knitr)
library(rlang)
library(stats)

#functions I need: for rounding, renaming columns in the data, my theme for graphs
mylog_trans <- function(base = exp(1), from = 0) {
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), 
            domain = c(base^from, Inf))
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

Rename_qPCR_colnames <- function(col) {
  # Fix cases like "A1: Sample 1"
  col <- gsub(": Sample [0-9].*$", "", col) 
  # Get indexes of columns that start with X (temp colums)
  temp_num <- which(grepl(pattern =  '^X', x =  col))
  # Get indexes of columns that contain the well info
  well_num <- which(!grepl(pattern =  '^X', x =  col))
  # change the X to well_Temp format
  col[temp_num] <- paste0(col[well_num], '_Temp')
  
  col[well_num] <- paste0(col[well_num], '_Fluo')
  return(col)
}

eli_th <- list(
  theme_minimal(), 
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11, colour = "black", face = "bold"), 
        axis.text.y = element_text(size = 11, colour = "black"), 
        axis.text.x = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size = 11), 
        plot.title = element_text(size = 18, hjust = 0, face = "bold"), 
        strip.text.x = element_text(size = 11, colour = "black"), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        legend.background = element_blank(), 
        legend.key = element_blank(), 
        panel.grid.minor = element_blank(),
  )
)

#establishing directories for file extraction

main_dir <- file.path("~/Desktop/MSc_Thesis")
assay_name <- "TriLineage_Diff1"
file_dir <- file.path(main_dir, assay_name)

#loading data for Plate 1
ct_file_1 <- file.path(file_dir, c("Plate1/20230426_Ct.txt"))

Plate1 <- read.delim(file = ct_file_1, header = TRUE, sep = '\t', dec = ",", 
                     stringsAsFactors = F, skip = 1, row.names = NULL) %>%
  dplyr::select(c("Pos", "Cp")) %>% 
  setNames(c("Well","CT"))

meta_data_path <- file.path(file_dir, "metadata.xlsx")

metadata1 <- read_excel(path = meta_data_path, sheet = 'Plate1' ) %>%
  setNames(c("Well", "Sample", "Genotype", "TimePoint", "Target"))

Plate1 <- left_join(metadata1, Plate1, by = "Well")
Plate1$Plate <- "Plate1"

#loading data for Plate 2
ct_file_2 <- file.path(file_dir, c("Plate2/20230427_Ct.txt"))

Plate2 <- read.delim(file = ct_file_2, header = TRUE, sep = '\t', dec = ",", 
                     stringsAsFactors = F, skip = 1, row.names = NULL) %>%
  dplyr::select(c("Pos", "Cp")) %>% 
  setNames(c("Well","CT"))

meta_data_path <- file.path(file_dir, "metadata.xlsx")

metadata2 <- read_excel(path = meta_data_path, sheet = 'Plate2' ) %>%
  setNames(c("Well", "Sample", "Genotype", "TimePoint", "Target"))

Plate2 <- left_join(metadata2, Plate2, by = "Well")
Plate2$Plate <- "Plate2"

#loading data for Plate 3
ct_file_3 <- file.path(file_dir, c("Plate3/20230428_Ct.txt"))

Plate3 <- read.delim(file = ct_file_3, header = TRUE, sep = '\t', dec = ",", 
                     stringsAsFactors = F, skip = 1, row.names = NULL) %>%
  dplyr::select(c("Pos", "Cp")) %>% 
  setNames(c("Well","CT"))

meta_data_path <- file.path(file_dir, "metadata.xlsx")

metadata3 <- read_excel(path = meta_data_path, sheet = 'Plate3' ) %>%
  setNames(c("Well", "Sample", "Genotype", "TimePoint", "Target"))

Plate3 <- left_join(metadata3, Plate3, by = "Well")
Plate3$Plate <- "Plate3"

#loading data for Plate 4

ct_file_4 <- file.path(file_dir, c("Plate4/20230502_Ct.txt"))

Plate4 <- read.delim(file = ct_file_4, header = TRUE, sep = '\t', dec = ",", 
                     stringsAsFactors = F, skip = 1, row.names = NULL) %>%
  dplyr::select(c("Pos", "Cp")) %>% 
  setNames(c("Well","CT"))

meta_data_path <- file.path(file_dir, "metadata.xlsx")

metadata4 <- read_excel(path = meta_data_path, sheet = 'Plate4' ) %>%
  setNames(c("Well", "Sample", "Genotype", "TimePoint", "Target"))

Plate4 <- left_join(metadata4, Plate4, by = "Well")
Plate4$Plate <- "Plate4"

#loading data for Plate 5

ct_file_5 <- file.path(file_dir, c("Plate5/20230503_Ct.txt"))

Plate5 <- read.delim(file = ct_file_5, header = TRUE, sep = '\t', dec = ",", 
                     stringsAsFactors = F, skip = 1, row.names = NULL) %>%
  dplyr::select(c("Pos", "Cp")) %>% 
  setNames(c("Well","CT"))

meta_data_path <- file.path(file_dir, "metadata.xlsx")

metadata5 <- read_excel(path = meta_data_path, sheet = 'Plate5' ) %>%
  setNames(c("Well", "Sample", "Genotype", "TimePoint", "Target"))

Plate5 <- left_join(metadata5, Plate5, by = "Well")
Plate5$Plate <- "Plate5"

#loading data for Plate 6

ct_file_6 <- file.path(file_dir, c("Plate6/2023_05_09_Ct.txt"))

Plate6 <- read.delim(file = ct_file_6, header = TRUE, sep = '\t', dec = ",", 
                     stringsAsFactors = F, skip = 1, row.names = NULL) %>%
  dplyr::select(c("Pos", "Cp")) %>% 
  setNames(c("Well","CT"))

meta_data_path <- file.path(file_dir, "metadata.xlsx")

metadata6 <- read_excel(path = meta_data_path, sheet = 'Plate6' ) %>%
  setNames(c("Well", "Sample", "Genotype", "TimePoint", "Target"))

Plate6 <- left_join(metadata6, Plate6, by = "Well")
Plate6$Plate <- "Plate6"

#loading data for Plate 7

ct_file_7 <- file.path(file_dir, c("Plate7/20230807_Ct.txt"))

Plate7 <- read.delim(file = ct_file_7, header = TRUE, sep = '\t', dec = ",", 
                     stringsAsFactors = F, skip = 1, row.names = NULL) %>%
  dplyr::select(c("Pos", "Cp")) %>% 
  setNames(c("Well","CT"))

meta_data_path <- file.path(file_dir, "metadata.xlsx")

metadata7 <- read_excel(path = meta_data_path, sheet = 'Plate7' ) %>%
  setNames(c("Well", "Sample", "Genotype", "TimePoint", "Target"))

Plate7 <- left_join(metadata7, Plate7, by = "Well")
Plate7$Plate <- "Plate7"

#making a dataframe qPCR with all Plates added


qPCR <- rbind(Plate1, Plate2, Plate3, Plate4, Plate5, Plate6, Plate7)


#adding replicates within the dataframe

qPCR <- qPCR %>%
  dplyr::arrange(Sample) %>%
  dplyr::arrange(Target) %>%
  dplyr::group_by(Sample, Target, Genotype, TimePoint, Plate) %>%
  dplyr::mutate(replicate = paste("r_", 1:2, sep = "")) |>
  ungroup()

#making sure that everything is numeric

qPCR$CT <- as.numeric(qPCR$CT)

### Melt Curve analysis - separate for each plate; needed to delete wells with multiple peaks/no peaks 

MeltCurves_plate7 <- file.path(file_dir, "Plate7/20230807_MeltCurves.txt")

melt_curves <- read_delim(file = MeltCurves_plate7 , delim = "\t",
                          escape_double = FALSE, trim_ws = TRUE, 
                          locale = locale(decimal_mark = ","), 
                          show_col_types = FALSE, 
                          name_repair = as_function(~ Rename_qPCR_colnames(col = .x) ) ) 

melt_curves |>
  pivot_longer(cols = everything(), names_sep = "_", names_to = c('well', 'feature'), values_to = "value")  |>
  arrange(well) |>
  pivot_wider(id_cols = well, names_from = feature, values_fn = list) |>
  unnest(cols = c(Temp, Fluo)) |>
  setNames(c('Well', 'Temperature', 'Fluorescence')) -> melt_curves

cols_to_keep <- c("Well", "Sample", "Target", "Genotype", "TimePoint")

qPCR |>
  subset(Plate == "Plate7") |>
  select(all_of(cols_to_keep)) |>
  left_join( y = melt_curves, by = join_by(Well) ) -> melt_profiles

melt_profiles |>
  group_by(Target) |>
  mutate(Peak = max(Fluorescence, na.rm = T ) ) |>
  mutate(T_at_peak = case_when(Peak == Fluorescence ~ Temperature)) |>
  ungroup() |>
  select(c("Sample", "Target", "Temperature", "T_at_peak") ) |>
  subset( !is.na( T_at_peak) ) -> peaks_max


filter(melt_profiles, Target == 'T' & Genotype == "NE") %>%
ggplot(melt_profiles) +
    aes(x = Temperature, y = Fluorescence, colour = Sample) +
    facet_grid(TimePoint+Genotype ~ Target, scales = "free_y" ) +
     geom_vline(data = peaks_max, aes(xintercept = T_at_peak),
                color = 'black', linetype = 'solid') +
     geom_point( size = 0.5, show.legend = F) +
     geom_path(show.legend = F) +
      scale_colour_brewer(type = "div", name = "Sample") +
     scale_x_continuous(n.breaks = 15) +
     guides(color = guide_legend(override.aes = list(size = 3) ) )+
     coord_cartesian(xlim = c(80, 90)) +
     labs(x = "Temperature (Celsius)", y = "Fluorescence A.U.") +
     theme_classic() + 
     eli_th 


# Filtering those wells for each plate and appending the 'filtered qPCR' to qPCR

Plate1_wells <- c('O17', 'G18', 'F6', 'C18', 'D9', 'G9', 'F10', 'E9', 'F9', 'A10')
Plate2_wells <- c('E23', 'E24', 'H4', 'P22', 'P21')
Plate4_wells <- c('E1', 'E2', 'H1', 'H2', 'C1', 'C2', 'L22', 'P22', 'N22')
Plate5_wells <- c('P21')
Plate6_wells <- c('J11', 'H11', 'O18', 'O24', 'O23', 'L19', 'L20', 'M20', 'M14')
Plate7_wells <- c('D18', 'D11', 'D17', 'D14', 'O14', 'D19', 'O16', 'O15', 'D4', 'D3')

filtered_qPCR1 <- subset(qPCR, Plate == 'Plate1' & Well != Plate1_wells)
filtered_qPCR2 <- subset(qPCR, Plate == 'Plate2' & Well != Plate2_wells)
filtered_qPCR3 <- subset(qPCR, Plate == 'Plate3' )
filtered_qPCR4 <- subset(qPCR, Plate == 'Plate4' & Well != Plate4_wells)
filtered_qPCR5 <- subset(qPCR, Plate == 'Plate5' & Well != Plate5_wells)
filtered_qPCR6 <- subset(qPCR, Plate == 'Plate6' & Well != Plate6_wells)
filtered_qPCR7 <- subset(qPCR, Plate == 'Plate7' & Well != Plate7_wells & Target != 'Hoxa1')

filtered_qPCR <- rbind(filtered_qPCR1, filtered_qPCR2, filtered_qPCR3, 
                       filtered_qPCR4, filtered_qPCR5, filtered_qPCR6,
                       filtered_qPCR7)

filtered_qPCR <- subset(filtered_qPCR, ! Target   %in% c('Tgfbi') )

qPCR <- filtered_qPCR


#ordering samples to appear in each graph 

samples_order <- c('ZWT', 'ZA1', 'ZA2', 'ZA3', 'ZO6', 'ZO7', 'ZO8', 'ZKO1', 'ZKO3', 'ZKO4')
qPCR$Sample <- factor(qPCR$Sample, levels = samples_order)

target_order <- c('Rex1', 'Pax6', 'Pou5f1' , 'Sox17', 'Foxa1','Gata4', 'T', 'Flk1', 
                  'Gsc','Nes', 'Nanog', 'Actb', 'Sox1','Essrb', 'Hoxa1', 'Twist', 'Tbp')

qPCR$Target <- factor(qPCR$Target, levels = target_order)

qPCR$TimePoint <- factor(qPCR$TimePoint,
                         levels = c("ESC", "EB", "CM", "NE", "EN"))

qPCR$Genotype <- factor(qPCR$Genotype,
                        levels = c("WT", "CSex4", "∆ex4", "KO"))

qPCR$Well <- factor(qPCR$Well, levels = unique(mixedsort(qPCR$Well)) )


qPCR <- subset(qPCR, Sample != 'ZKO1') #this clone was not a KO in the end - we are excluding it from our analysis  

#plotting the Ct values for each replicate 
ggplot(qPCR, aes(x = Sample, y = CT, group = replicate, shape = factor(TimePoint), fill = Target)) +
  facet_wrap(~ Target, scales = "free", ncol = 5, nrow = 5) +
  geom_point(size = 1.5, alpha = 0.65, 
             position = position_dodge(width = 0.5), 
             show.legend = TRUE) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "TimePoint") +
  guides(fill = "none") +
  coord_cartesian(clip = 'off', ylim = c(15, 35)) +
  labs(title = "Ct values in each well", x = "", y = "Cycle Number") +
  eli_th

write_delim(x = qPCR, file = file.path(file_dir, 'Combined_raw_ct_values.tab'), 
            delim = '\t', append = F, quote = 'none')

#seeing the plate effect 

qPCR |>
  group_by(Plate) |>
  mutate(Rows = str_split_fixed(Well, pattern = "", n = 2)[,1],
         Cols = str_split_fixed(Well, pattern = "", n = 2)[,2],
         .after = 'Well') |>
  mutate(Cols = as.integer(Cols)) |>
  droplevels() |>
  ggplot() +
  aes(x = Cols, y = reorder(Rows, desc(Rows) ) , fill = CT) +
  facet_wrap(~Plate) +
  geom_tile() +
  labs(x = "", y = "") +
  guides(fill = guide_colourbar(barwidth = unit(x = 12, 'cm'))) +
  scale_fill_viridis(name = "CT value", breaks = seq(15, 45, 2), direction = -1) +
  scale_x_continuous(n.breaks = 24) +
  coord_cartesian(expand = F) + 
  eli_th +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

df_reps_sum <- qPCR %>%
  group_by(Sample, Target, TimePoint) %>%
  dplyr::summarise(CT_mean = mean(CT, na.rm = TRUE), 
                   CT_StDev = sd(CT, na.rm = TRUE) ) %>%
  ungroup()


qPCR <- left_join(qPCR, df_reps_sum, 
                  by = c("Sample", "Target", "TimePoint"))
#seeing the StDev of each sample 

ggplot(qPCR) +
  aes(x = Sample, y = CT_StDev, shape = factor(TimePoint)) +
  facet_wrap( ~ Target, scale = "fixed", nrow = 1) +
  geom_point(aes(fill = CT_StDev < 0.65), size = 2, show.legend = T,
             position = position_dodge(width = 0.5) ) +
  geom_hline(yintercept = 1, color = "gray16", linetype = 6 ) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "Time Point") +
  guides(fill = "none") +
  scale_fill_manual(values = c('aquamarine4', 'coral2'), guide = F) +
  scale_y_continuous(trans = mylog_trans(base = 10, from = -2.0), 
                     breaks = c(0, 0.015, 0.03, 0.062, 0.125, 0.25, 0.5, 1, 2), 
                     limits = c(NA, 2), 
                     oob = squish) + 
  labs(y = "Standard deviation", x = "", caption = "y-axis is in log10 scale") +
  eli_th

###Analysis part - comparative method using the formula ∆∆CT (GOI CT - Housekeeping CT) - (Control ∆CT), then negative of this is the exponent of base 2
###Housekeeping gene is Tbp 

HKG <- subset(qPCR, Target %in% c("Tbp")) %>% #  c("Tbp", "Actb")
  group_by(Sample, Target, TimePoint) %>%
  dplyr::summarize(CT_mean_norm = mean(CT, na.rm = TRUE), 
                   CT_StDev_norm = sd(CT, na.rm = TRUE)) %>%
  ungroup()

HKG$Target <- NULL

ggplot(HKG) +
  aes(x = Sample, y = CT_StDev_norm, shape = factor(TimePoint)) +
  geom_point(aes(fill = CT_StDev_norm < 0.65), size = 2, show.legend = T,
             position = position_dodge(width = 0.5) ) +
  geom_hline(yintercept = 1, color = "gray16", linetype = 6 ) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "Time Point") +
  guides(fill = "none") +
  scale_fill_manual(values = c('aquamarine4', 'coral2'), guide = F) +
  scale_y_continuous(trans = mylog_trans(base = 10, from = -2.0), 
                     breaks = c(0, 0.015, 0.03, 0.062, 0.125, 0.25, 0.5, 1, 2), 
                     limits = c(NA, 2), 
                     oob = squish) + 
  labs(y = "Standard deviation", x = "", caption = "y-axis is in log10 scale") +
  eli_th

qPCR <- left_join(qPCR, HKG, by = c("Sample", "TimePoint")) %>%
  mutate(Delta_CT = CT - CT_mean_norm,
         Delta_CT_StDev = sqrt((CT_StDev ^2) + (CT_StDev_norm ^2)))


#∆∆CT to WT 
#First create the reference ∆CT, and add back to the datafram and calculate the ∆∆CT
ref_Delta_CT <- subset(qPCR, Genotype == "WT") %>%
  dplyr::select(Target, replicate, TimePoint, Delta_CT) %>%
  group_by(Target, TimePoint) %>%
  dplyr::summarize(ref_Delta_CT = mean(Delta_CT, na.rm = TRUE))

qPCR_WT <- left_join(qPCR, ref_Delta_CT, by = c("Target", "TimePoint") ) %>%
  dplyr::mutate(Delta_Delta_CT = Delta_CT - ref_Delta_CT)

#fold change calculation
qPCR_WT <-
  group_by(qPCR_WT, Sample, Target, TimePoint) %>%
  dplyr::summarize(Delta_Delta_CT_mean = mean(Delta_Delta_CT, na.rm = TRUE), 
                   Delta_Delta_CT_std = sqrt(var(Delta_Delta_CT, na.rm = TRUE))) %>% 
  ungroup() %>% 
  left_join(qPCR_WT, ., by = c("Sample", "Target", "TimePoint") ) %>%
  dplyr::mutate(Comp_Meth = 2^-(Delta_Delta_CT_mean),
                Comp_Meth_SD.plu = 2^ - ( Delta_Delta_CT_mean - Delta_CT_StDev),
                Comp_Meth_SD.min = 2^ - ( Delta_Delta_CT_mean + Delta_CT_StDev) )

# sample ~1% of the dataframe to print
sample_frac(round_df(qPCR_WT, 2), size = 0.4)

#Plotting the StDev of the ∆∆CT

ggplot(qPCR_WT) +
  aes(x = Sample, 
      y = 2^-Delta_Delta_CT_mean, 
      ymin = Comp_Meth_SD.min,  
      ymax = Comp_Meth_SD.plu, 
      fill = Target)  +
  facet_grid(Target ~ TimePoint, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = 6, color = "grey84") +
  geom_col(position = position_dodge(), width = 0.5,
           alpha = 0.4, colour = "black", show.legend = F) +
  geom_errorbar(stat = "identity", width = 0.25, 
                position = position_dodge(width = 0.55), color = "black", 
                alpha = 0.25) +
  geom_point(inherit.aes = F, show.legend = F,
             aes(x = Sample, y = 2^-Delta_Delta_CT, 
                 fill = Target, group = replicate), 
             shape = 21, col = "black", size = 1.5, alpha = 0.85, 
             position = position_dodge(width = 0.35) ) +
  guides(fill = guide_legend(nrow = 1) ) +
  scale_y_continuous(expand = expansion(add = c(0.01, 0.5)),
                     trans = mylog_trans(base = 2, from = -8.5),
                     labels = scales::number_format(accuracy = 0.01),
                     n.breaks = 6) +
  labs(x = "Sample", 
       y = "Fold Change\n(normalised to Tbp and WT\at each developmental stage)",
       title = "", caption = "Y-axis on log2 scale") +
  eli_th + theme(panel.grid.major.x = element_blank())

#Plotting as a heatmap to WT cells

qPCR_WT |>
  select(Well, Sample, Genotype, TimePoint, Target,  Plate, Delta_Delta_CT_mean) |>
  group_by(TimePoint, Target, Plate, Genotype) |>
  mutate(Genotype_Mean_Delta_Delta_CT = mean(Delta_Delta_CT_mean, na.rm = T),
         Genotype_Mean_Delta_Delta_CT = round(Genotype_Mean_Delta_Delta_CT, 1) ) |>
  select(Genotype, TimePoint, Target,  Plate, Genotype_Mean_Delta_Delta_CT) |>
  unique() -> qPCR_WT_Mean_by_Genotype

#defining a new target order to have pluripotency on top, then lineage specific 

new_target_order <- c("Nanog", "Rex1", "Pax6", "Foxa1", "Gata4",
                      "Flk1", "Actb", "Sox1", "Hoxa1", "Twist")
  
qPCR_WT_Mean_by_Genotype |>
  subset(Target %in% c("Nanog", "Rex1", "Pax6", "Foxa1", "Gata4",
                       "Flk1", "Actb", "Sox1", "Hoxa1", "Twist")) |> 
  mutate(Target = factor(Target, levels = new_target_order)) |>
  ggplot() +
      aes(x = TimePoint, 
          y = Genotype, 
          fill = 2^-Genotype_Mean_Delta_Delta_CT, 
          title("Pluripotency Gene Expression Change (in comparison to WT)")) +
  facet_wrap( ~ Target, ncol = 4)+
  geom_tile(colour = 'black', linewidth = 0.1) +
  geom_text(aes(label = round(2^-Genotype_Mean_Delta_Delta_CT, 1) ) ) +
  scale_fill_gradient2(low = 'dodgerblue', mid = 'white', midpoint = 1, high = 'firebrick3', 
                       na.value = 'gray', limits = c(0, 5), oob = squish ) +
  theme_bw(base_size = 12) +
  labs(x = "", 
       y = "", caption = "ESC - Embryonic Stem Cells; EB - Embryoid Body; CM - Cardio-Mesoderm; EN - Endoderm; NE - Neuro-Ectoderm", 
       fill = "log2 fold change") +
  theme(legend.position = 'bottom')+
  eli_th


#Normalisation to ESCs, and added back into qPCr dataset; calculation of ∆∆CT (makes a new column in qPCR)

ref_Delta_CT_ESC <- subset(qPCR, Genotype == "WT" & TimePoint == 'ESC') %>%
  select(Target, replicate, Delta_CT) %>%
  group_by(Target) %>%
  dplyr::summarize(ref_Delta_CT = mean(Delta_CT, na.rm = TRUE) )

qPCR_ESC <- left_join(qPCR, ref_Delta_CT_ESC, by = c("Target") ) %>%
  dplyr::mutate(Delta_Delta_CT = Delta_CT - ref_Delta_CT)

# now fold change, first calculate the ∆∆CT between clones of the same genotype, and then use the formula for WT comparison
qPCR_ESC <-
  group_by(qPCR_ESC, Sample, Target, TimePoint) %>%
  dplyr::summarize(Delta_Delta_CT_mean = mean(Delta_Delta_CT, na.rm = T), 
                   Delta_Delta_CT_std = sqrt(var(Delta_Delta_CT, na.rm = T))) %>% 
  ungroup() %>% 
  left_join(qPCR_ESC, ., by = c("Sample", "Target", "TimePoint") ) %>%
  dplyr::mutate(Comp_Meth = 2^-(Delta_Delta_CT_mean),
                Comp_Meth_SD.plu = 2^ - ( Delta_Delta_CT_mean - Delta_CT_StDev),
                Comp_Meth_SD.min = 2^ - ( Delta_Delta_CT_mean + Delta_CT_StDev))


# sample ~1% of the dataframe to print
sample_frac(round_df(qPCR_ESC, 2), size = 0.4)

#Plotting as a heatmap to ESC

qPCR_ESC |>
  select(Well, Sample, Genotype, TimePoint, Target,  Plate, Delta_Delta_CT_mean) |>
  group_by(TimePoint, Target, Plate, Genotype) |>
  mutate(Genotype_Mean_Delta_Delta_CT = mean(Delta_Delta_CT_mean, na.rm = T),
         Genotype_Mean_Delta_Delta_CT = round(Genotype_Mean_Delta_Delta_CT, 1) ) |>
  select(Genotype, TimePoint, Target,  Plate, Genotype_Mean_Delta_Delta_CT) |>
  unique() -> qPCR_ESC_Mean_by_Genotype


qPCR_ESC_Mean_by_Genotype |>
  subset(Target %in% c("Nanog", "Rex1", "Pax6", "Foxa1", "Gata4",
                       "Flk1", "Actb", "Sox1", "Hoxa1", "Twist")) |> 
  mutate(Target = factor(Target, levels = new_target_order)) |>
  filter (Target != 'Tbp') |>
  ggplot()+
      aes(x = TimePoint, 
          y = Genotype, 
          fill = 2^-Genotype_Mean_Delta_Delta_CT)+
          #title = "Pluripotency Gene Expression Change (in comparison to ESC)") +
  facet_wrap( ~ Target, ncol = 4)+
  geom_tile(colour = 'black', linewidth = 0.1) +
  geom_text(aes(label = round(2^-Genotype_Mean_Delta_Delta_CT, 1) ) ) +
  scale_fill_gradient2(low = 'dodgerblue', mid = 'white', midpoint = 1, high = 'firebrick3', 
                       na.value = 'gray', limits = c(0, 15), oob = squish ) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom')+
  labs(x = "", 
       y = "", caption = "ESC - Embryonic Stem Cells; EB - Embryoid Body; CM - Cardio-Mesoderm; EN - Endoderm; NE - Neuro-Ectoderm", 
       fill = "log2 fold change") +
  theme(legend.position = 'bottom')+
  eli_th
