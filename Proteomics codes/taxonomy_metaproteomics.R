# Analysis of metaproteomics dataset from human colon 
# Written by Johan Sebastián Sáenz https://github.com/SebasSaenz
# Univeristy of Hohenheim, Department of Functional Microbiology of Livestock
# https://livestock-functional-microbiology.uni-hohenheim.de/en/english

# This code analyze the abundance and structure of peptides (lfq) assigned to
# different bacterial taxonomic levels. Data set was obtained after processing
# raw metaproteomic files with MetaLab.

################################################################################

### Load libraries  using pacman ###
# If not installed libraries will be install
pacman::p_load(tidyverse,
               ggtext,
               purrr,
               broom,
               future.apply,
               patchwork,
               vegan)

### set working directory ###
setwd("~/NextCloud/Rusitec_data")

### load files ###
df_taxonomy <- read.table("Raw_data/BuiltIn.taxa.refine.csv", # taxa abundance
                     sep = ",",
                     header = TRUE,
                     check.names = FALSE)

df_metadata <- read.csv("Raw_data/hHSILnew_mod.csv", #deriv cohort metadata
                    sep = "\t") %>%
  filter(deriv == 1)  %>%
  filter(bHSILnew != "") %>%
  mutate(record_id = str_replace(record_id, "-", "_"),
        bHSILnew = str_replace(bHSILnew, "bHSIL", "HSIL")) %>% 
  select(record_id, bHSILnew)

df_meta_all <- read.csv("Raw_data/all_metadata.txt", #all samples metadata
                        sep = "\t") %>%
  select(record_id, cAINcat) %>%
  mutate(record_id = str_replace(record_id, "-", "_"))

### Ordination NMDS ###
df_taxonomy$row_num <- seq.int(nrow(df_taxonomy)) #make each row unique

nmds_df <- df_taxonomy %>%  # make a long data frame
  pivot_longer(cols = c(11:290), 
               names_to = "record_id", 
               values_to = "lfq") %>%
  select(row_num, everything()) %>%
  inner_join(df_metadata, by = "record_id") %>%
  pivot_wider(id_cols = c(1:11), names_from = "record_id", values_from = "lfq")
  
nmds_df <- nmds_df[c(12:158)] %>% #calculate relative abundance and create a matrix 
  apply(.,2, function(x){(x/sum(x))*100}) %>%
  t() 

nmds1 <- metaMDS(nmds_df,  #perform nmds
                 distance = "bray", 
                 try = 20, 
                 trymax = 100, 
                 maxit = 1000, 
                 k = 3)


nmds_best <- metaMDS(nmds_df, distance = "bray", #find best nmds
                     try = 20, 
                     trymax = 100, 
                     maxit = 1000, 
                     k = 3,
                     previous.best = nmds1)

# fast ploting and qc checking
nmds_best$stress  
stressplot(nmds_best)
plot(nmds_best)


#plotting with ggplot
data.scores <- as.data.frame(scores(nmds_best, display=c("sites")))

#Addd metadata to dataframe
data.scores$record_id <- as.character(row.names(data.scores)) 

#joing metadata nmds scores
data_nmds <- left_join(data.scores, df_metadata, by = "record_id")

data_nmds$bHSILnew <- factor(data_nmds$bHSILnew,
                             levels = c("No HSIL", "HSIL"))
#make plot bHSIL
nmds <- data_nmds %>%
  ggplot() +
  stat_ellipse(linetype = 2,
               aes(x = NMDS1,
                   y = NMDS2,
                   colour = bHSILnew),
               level = 0.95) +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 colour =bHSILnew),
             size = 2
             #alpha = 0.7
             ) +
  scale_color_manual(values = c('#66C2A5','#FC8D62')) +
  theme_classic() +
  theme(panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "bottom") #remove legend title

#save figure
ggsave(nmds, filename = "Figures/nmds_bHSIL.png", width = 4, height = 4)

# ordination with cAINcat

#join the whole metadata witn nmds scores
caincat_meta <- inner_join(df_metadata, df_meta_all, by = "record_id")

data_nmds_caincat <- left_join(data.scores, caincat_meta, by = "record_id")

data_nmds_caincat$cAINcat <- factor(data_nmds_caincat$cAINcat,
                             levels = c("Normal", "LSIL", "HSIL"))
# Make plot
nmds <- data_nmds_caincat %>%
  ggplot() +
  stat_ellipse(linetype = 2,
               aes(x = NMDS1, 
                   y = NMDS2,
                   colour = cAINcat),
                  level = 0.95) +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 colour =cAINcat),
             size = 2
             #alpha = 0.7
             ) +
  scale_color_manual(values = c('#8DA0CB','#66C2A5','#FC8D62')) +
  theme_classic() +
  theme(panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "bottom") #remove legend title


ggsave(nmds, filename = "Figures/nmds_caincat.png", width = 4, height = 4)

#### statistics ordination (PERMANOVA) ###

# Arrange data to perfrom PERMANOVA test
nmds_df_stats <- df_taxonomy %>%
  pivot_longer(cols = c(11:290), 
               names_to = "record_id", 
               values_to = "lfq") %>%
  select(row_num, everything()) %>%
  inner_join(caincat_meta, by = "record_id") %>%
  pivot_wider(id_cols = c(1:11), names_from = "record_id", values_from = "lfq")

# Make a metadata for the test
nmds_meta <- df_taxonomy %>%
  pivot_longer(cols = c(11:290), names_to = "record_id", values_to = "lfq") %>%
  select(row_num, everything()) %>%
  inner_join(caincat_meta, by = "record_id") %>%
  select(record_id, bHSILnew, cAINcat) %>%
  unique()

# Calculate relative abundance for this new df
stats_df <- nmds_df_stats[c(12:158)] %>%
  apply(.,2, function(x){(x/sum(x))*100}) %>%
  t() 

# Calculate distance
tax.dist <- vegdist(stats_df, method = "bray") 

#Run PERMANOVA
adonis2(tax.dist ~ cAINcat, 
        data = nmds_meta, 
        permutations = 999)

########################## Barplots #########################################

# You can mofdify to plot bHSIL or cAINcat

pool <- 0.1 # threshold to pool taxonmic level
  
  rel_abundance <<- df_taxonomy %>%
    select(-c(1:4)) %>%
    pivot_longer(cols = c(7:286), names_to = "record_id", values_to = "lfq") %>%
      select(record_id, Phylum, lfq) %>% #select tax level
    mutate(Phylum = str_replace(Phylum, "(\\S*)$", "*\\1*")) %>% # for italics
    mutate(Phylum = str_replace(Phylum, "^\\**$", "")) %>%
    mutate(Phylum = str_replace(Phylum, "^$", "Unclassified")) %>% # rename empty 
    group_by(Phylum, record_id) %>%
    summarize(lfq = sum(lfq), .groups="drop") %>%
    group_by(record_id) %>%
    mutate(rel_abund = 100*(lfq/sum(lfq))) %>%
    ungroup() %>%
    inner_join(caincat_meta, by = "record_id") %>% # join metadata
    select(Phylum, record_id, rel_abund, bHSILnew) %>% # change bHSILnew/cAINcat
    filter(Phylum != "")
  
  
  Phylum_pool <<- rel_abundance %>% # pool tax level absed on treshold
    group_by(Phylum) %>%
    summarize(pooled = mean(rel_abund) < pool,
              mean = mean(rel_abund),
              .groups ="drop")
  
  rel_abundance_pool <<- rel_abundance %>% #mean relative abundance
    inner_join(Phylum_pool, by = "Phylum") %>%
    mutate(Phylum=if_else(pooled, "Other", Phylum)) %>%
    group_by(Phylum, bHSILnew) %>%
    summarise(mean_rel_abund = mean(rel_abund), .groups = "drop") 

# set colors
col_list <- c('#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69',
              '#fccde5','#d9d9d9', '#8dd3c7','#bc80bd','#ccebc5','#ffed6f',
              "black")

# plot
 barplot_tax <-  rel_abundance_pool %>%
    filter(Phylum != "Unclassified")%>% # remove Unclassified
    ggplot(aes(x = bHSILnew,
               y = mean_rel_abund,
               fill = Phylum)) +
    geom_bar(stat = "identity", 
             width = 0.6) +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,70)) +
    scale_fill_manual(values = col_list) +
    labs(x = "",
         y = "Relative abundance (%)") +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.text = element_markdown(),
          legend.key.size = unit(11, "pt"),
          axis.text.x = element_text(size = 12))

ggsave(barplot_tax, 
       filename = "Figures/barplot_phylum_bHSILnew.png", 
       width = 5, 
       height = 4)

#### Kruskal-wallis test for relative abundance ###

significant <- rel_abundance %>%
    nest(data = -Phylum) %>%
    mutate(test = map(.x = data, ~kruskal.test(rel_abund~bHSILnew, data = .x) %>%
                        tidy())) %>%
    unnest(test) %>%
    mutate(p.adjusted = p.adjust(p.value, method = "BH")) %>%
    filter(p.adjusted < 0.05) %>% # filter by adjusted p-value
    select(Phylum, p.adjusted)


sort_samples <- rel_abundance %>% # sort samples
  group_by(Phylum) %>%
  summarise(mean_rel = mean(rel_abund)) %>%
  arrange(desc(mean_rel)) %>%
  pull(Phylum)

# make factors
rel_abundance$bHSILnew <- factor(rel_abundance$bHSILnew,
                                    levels = c("HSIL", "No HSIL"))


# lollipop plot
lollipop <- rel_abundance %>%
    mutate(Phylum = factor(Phylum, levels = sort_samples),
           bHSILnew = factor(bHSILnew, levels = c("HSIL", "No HSIL")))%>%
    filter(Phylum != "Unclassified", # filter unwanted groups
           Phylum != "*Nematoda*",
           Phylum != "*Chordata*") %>%
    ggplot(aes(x = rel_abund+0.0001,
               y = Phylum,
               color = bHSILnew,
               fill = bHSILnew)) +
    geom_jitter(size = 1,
                position = position_jitterdodge(dodge.width = 0.8,
                                                jitter.width = 0.3),
                shape = 21) +
    stat_summary(size = 0.3,
                 fun.data = median_hilow, 
                 fun.args = list(conf.int = 0.50),
                 geom = "pointrange",
                 position = position_dodge(width = 0.8),
                 color = "black",
                 show.legend = FALSE) +
    labs(x = "", 
         y ="Relative abundance (%)") +
    theme_classic()  +
    theme(axis.text.y  = element_markdown(),
          legend.title = element_blank(),
          legend.position = "bottom") +
    scale_x_log10(labels = function(x) format(x, scientific = FALSE), #rename log scale
                                              breaks = c(0, 0.0006, 0.006, 0.06, 0.6, 6, 60))+
    scale_fill_manual(NULL,
                      breaks = c("HSIL", "No HSIL"),
                      values = c('#FC8D62','#66C2A5')) +
    scale_color_manual(NULL,
                       breaks = c("HSIL", "No HSIL"),
                       values = c('#FC8D62','#66C2A5'))

ggsave(filename = "Figures/dotplot_phylum_bhsil.png", width = 4, height = 5)



#test specific groups

filtered_tax <- rel_abundance %>%
  filter(Phylum == "*Bacteroidetes*") #filter group

kruskal.test(rel_abund~bHSILnew, data = filtered_tax) #test

# plot 

filtered_tax %>% 
  ggplot(aes(y = rel_abund,
             x = bHSILnew,
             color = bHSILnew)) + 
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.6)) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.95),
               geom = "pointrange",
               position = position_dodge(width = 0.8),
               size=0.5,
               color="black") +
  labs(y = "Relative abundance (%)",
       x = "") +
  theme_classic()
  
  



