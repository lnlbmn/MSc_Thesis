library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstatix)

eli_th <- list(
  theme_bw(), 
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
        strip.background = element_blank()
  )
)

oneDrive_dir <- file.path('~/OneDrive - CRG - Centre de Regulacio Genomica/Suz12_AS_project')
resazurin_dir <- file.path(oneDrive_dir, '04_ESCs_Differentiations/09_Resazurin_Differentiation_noLIF')

raw_data <- read_excel(path = file.path(resazurin_dir, "Data_23_04_19.xlsx"),
                       sheet = "500_cells")

# to reshape to long format
raw_data |> 
  pivot_longer(cols = !c(CellNumber, ReadingTime, WaveLength, TimePoint), 
               names_to = 'CellLine', values_to = 'Absorbance') -> long_data 

# the molar extinction coefficient of the oxidized resazurin at 600nm: 117216
Eoxi600 <- 117216

# the molar extinction coefficient of the oxidized resazurin at 570nm: 80586
Eoxi570 <- 80586

# the molar extinction coefficient of the reduced resazurin at 600nm: 
Ered600 <- 14652

# the molar extinction coefficient of the reduced resazurin at 600nm: 
Ered570 <- 155677

# normalisation of the denominator
long_data |>
  group_by(CellNumber, ReadingTime, TimePoint) |>
  reframe(Denominator = (Ered570*Absorbance[WaveLength == 600 & CellLine == 'empty']) - (Ered600*Absorbance[WaveLength == 570 & CellLine == 'empty']) ) |>
  group_by(CellNumber, ReadingTime, TimePoint) |>
  summarise(AverageCntrl = mean(Denominator) ) -> normaliser
merged_data <- left_join(long_data, y = normaliser, by = join_by(CellNumber, ReadingTime, TimePoint))

merged_data |> 
  subset(WaveLength %in% c(570, 600)) |>
  group_by(ReadingTime, TimePoint) |>
  mutate(percent = ( ( (Eoxi600*Absorbance[WaveLength == 570]) - (Eoxi570*Absorbance[WaveLength == 600]) ) / AverageCntrl )*100 ) -> norm_data

annotation_data <- data.frame(CellLine = c("WT", "ZA1", "ZA2", 
                                           "ZA3I", "ZA3II",
                                           "ZO6", "ZO7", "ZO8",
                                           "ZKO1", "ZKO3", "ZKO4", "empty"), 
                              Genotype = c('WT', 'WT', 'WT', 'CSex4', 'CSex4', 
                                           '∆ex4', '∆ex4', '∆ex4', 
                                           'KO', 'KO', 'KO', 'empty'))
annotation_data$Genotype <- factor(annotation_data$Genotype, levels = unique(annotation_data$Genotype))

norm_data <- left_join(norm_data, annotation_data, by = join_by(CellLine))

ggplot(norm_data) +
  aes(x = Genotype, y = percent, fill = Genotype) +
  facet_grid(CellNumber~ReadingTime+TimePoint, scales = 'free_y') +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 21, show.legend = F) +
  # coord_cartesian(ylim = c(0, 100)) +
  labs(y = '% Reduction of resazurin') +
  theme(axis.text = element_text(colour = 'black'),
        legend.position = 'none')
#ggsave(filename = 'Resazurin_percentage.pdf', device = 'pdf', path = dwnfdr, 
       #units = 'cm', width = 18, height = 5)

unique(norm_data$ReadingTime)
norm_data |>
  ungroup() |>
  group_by(CellLine, ReadingTime) |>
  mutate(Rel_percent = percent / percent[TimePoint == 0] ) |>
  # mutate(Rel_percent = percent - percent[Genotype == 'WT'] ) |>
  subset(Genotype != 'empty') |>
  #subset(CellNumber == 500) |>
  subset(ReadingTime != 180) -> tidy_norm_data


library(ggsignif)


ggplot(tidy_norm_data) +
  aes(x = Genotype, y = Rel_percent, fill = Genotype) +
  facet_grid(CellNumber ~ TimePoint) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape = factor(ReadingTime)), show.legend = F) +
  scale_shape_manual(values = c("80" = 21, "120" = 22)) + #21 - circle, 22 - square 
  scale_fill_manual(values = c('WT' = "#377AA3", '∆ex4' = "#F7CB48",
                               'CSex4' = "#B65120",'KO' = "#728189", 
                               'KO+L' = "#7B67AB", 'KO+S' = "#F07E19") )  +
  geom_hline(yintercept = 1) +
  labs(y = '% Reduction of resazurin relative to day zero') +
  eli_th +
  theme(axis.text = element_text(colour = 'black'),
        legend.position = 'none')

#tidy_norm_data$Genotype <- factor(tidy_norm_data$Genotype, levels = c("WT", "CSex4", "∆ex4", "KO"))
#unique(tidy_norm_data$Genotype)
#dput(head(tidy_norm_data,10))


tidy_norm_data |>
  filter(TimePoint == 4) |> 
  filter(WaveLength == 600) |>
  filter(ReadingTime == 120) |> 
  arrange(CellLine) |>
  droplevels() |>
  ungroup() |>
  wilcox_test(percent ~ Genotype, ref.group = 'WT', 
              paired = F, exact = T, alternative = "two.sided",
              conf.level = 0.95, detailed = T) |>
  adjust_pvalue(method = "hochberg") |>
  add_significance("p.adj")

library(ggsignif)


ggplot(tidy_norm_data) +
  aes(x = Genotype, y = Rel_percent, fill = Genotype) +
  facet_grid(CellNumber ~ TimePoint) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape = factor(ReadingTime)), show.legend = F) +
  scale_shape_manual(values = c("80" = 21, "120" = 22)) + #21 - circle, 22 - square 
  scale_fill_manual(values = c('WT' = "#377AA3", '∆ex4' = "#F7CB48",
                               'CSex4' = "#B65120",'KO' = "#728189", 
                               'KO+L' = "#7B67AB", 'KO+S' = "#F07E19") )  +
  geom_hline(yintercept = 1) +
  labs(y = '% Reduction of resazurin relative to day zero') +
  eli_th +
  geom_signif(
    comparisons = list(c("WT", "CSex4"), c("WT", "∆ex4"), c("WT", "KO")) +
      map_signif_level = TRUE, size = 6) +
  theme(axis.text = element_text(colour = 'black'),
        legend.position = 'none')
