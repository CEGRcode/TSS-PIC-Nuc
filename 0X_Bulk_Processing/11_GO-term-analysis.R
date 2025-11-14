
## in R
library(rGREAT)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(forcats)
getwd()
setwd("/Path/to/Title/Library/E11")
list.files()

bed_YRWS <- read.table("./TSS_+1Nuc_YRWS.bed")
gr_YRWS <- GRanges(seqnames = bed_YRWS[, 1], ranges = IRanges(bed_YRWS[, 2], bed_YRWS[, 3]))
bed_bg <- read.table("./TSS_+1Nuc_all.bed")
gr_bg <- GRanges(seqnames = bed_bg[, 1], ranges = IRanges(bed_bg[, 2], bed_bg[, 3]))

res_YRWS <- great(gr_YRWS, "GO:BP", "Gencode_v38", background = gr_bg)
tb3 <- getEnrichmentTable(res_YRWS) %>% mutate(source = "YRWS")
write.csv(tb3, "tb3_results.csv", row.names = FALSE)

bed_lesscode <- read.table("./TSS_+1Nuc_lessDNAencode.bed")
gr_lesscode <- GRanges(seqnames = bed_lesscode[, 1], ranges = IRanges(bed_lesscode[, 2], bed_lesscode[, 3]))

res_lesscode <- great(gr_lesscode, "GO:BP", "Gencode_v38", background = gr_bg)
tb9 <- getEnrichmentTable(res_lesscode) %>% mutate(source = "lesscode")
write.csv(tb9, "tb9_results.csv", row.names = FALSE)


bed_sameYR_lowWS <- read.table("./TSS_+1Nuc_sameYR_lowWS.bed")
gr_sameYR_lowWS <- GRanges(seqnames = bed_sameYR_lowWS[, 1], ranges = IRanges(bed_sameYR_lowWS[, 2], bed_sameYR_lowWS[, 3]))

res_sameYR_lowWS <- great(gr_sameYR_lowWS, "GO:BP", "Gencode_v38", background = gr_bg)
tb4 <- getEnrichmentTable(res_sameYR_lowWS) %>% mutate(source = "sameYR_lowWS")
write.csv(tb4, "tb4_results.csv", row.names = FALSE)

bed_lowYR_sameWS <- read.table("./TSS_+1Nuc_lowYR_sameWS.bed")
gr_lowYR_sameWS <- GRanges(seqnames = bed_lowYR_sameWS[, 1], ranges = IRanges(bed_lowYR_sameWS[, 2], bed_lowYR_sameWS[, 3]))

res_lowYR_sameWS <- great(gr_lowYR_sameWS, "GO:BP", "Gencode_v38", background = gr_bg)
tb5 <- getEnrichmentTable(res_lowYR_sameWS) %>% mutate(source = "lowYR_sameWS")
write.csv(tb5, "tb5_results.csv", row.names = FALSE)

bed_sameYR_antiWS <- read.table("./TSS_+1Nuc_sameYR_antiWS.bed")
gr_sameYR_antiWS <- GRanges(seqnames = bed_sameYR_antiWS[, 1], ranges = IRanges(bed_sameYR_antiWS[, 2], bed_sameYR_antiWS[, 3]))

res_sameYR_antiWS <- great(gr_sameYR_antiWS, "GO:BP", "Gencode_v38", background = gr_bg)
tb6 <- getEnrichmentTable(res_sameYR_antiWS) %>% mutate(source = "sameYR_antiWS")
write.csv(tb6, "tb6_results.csv", row.names = FALSE)

bed_antiYR_sameWS <- read.table("./TSS_+1Nuc_antiYR_sameWS.bed")
gr_antiYR_sameWS <- GRanges(seqnames = bed_antiYR_sameWS[, 1], ranges = IRanges(bed_antiYR_sameWS[, 2], bed_antiYR_sameWS[, 3]))

res_antiYR_sameWS <- great(gr_antiYR_sameWS, "GO:BP", "Gencode_v38", background = gr_bg)
tb7 <- getEnrichmentTable(res_antiYR_sameWS) %>% mutate(source = "antiYR_sameWS")
write.csv(tb7, "tb7_results.csv", row.names = FALSE)

bed_antiYR_antiWS <- read.table("./TSS_+1Nuc_antiYR_antiWS.bed")
gr_antiYR_antiWS <- GRanges(seqnames = bed_antiYR_antiWS[, 1], ranges = IRanges(bed_antiYR_antiWS[, 2], bed_antiYR_antiWS[, 3]))

res_antiYR_antiWS <- great(gr_antiYR_antiWS, "GO:BP", "Gencode_v38", background = gr_bg)
tb8 <- getEnrichmentTable(res_antiYR_antiWS) %>% mutate(source = "antiYR_antiWS")
write.csv(tb8, "tb8_results.csv", row.names = FALSE)


tb9 <- read.csv("tb9_results.csv")
tb8 <- read.csv("tb8_results.csv")
tb7 <- read.csv("tb7_results.csv")
tb6 <- read.csv("tb6_results.csv")
tb5 <- read.csv("tb5_results.csv")
tb4 <- read.csv("tb4_results.csv")
tb3 <- read.csv("tb3_results.csv")


functions_of_interest <- c(
  "regulation of protein transport",
  "natural killer cell proliferation",
  "positive regulation of insulin receptor signaling pathway",
  "regulation of establishment of protein localization",
  "endoplasmic reticulum mannose trimming",
  "calcium-mediated signaling",
  "mRNA metabolic process",
  "extrinsic apoptotic signaling pathway via death domain receptors",
  "extracellular matrix assembly",
  "negative regulation of cell junction assembly",
  "vitamin transport",
  "response to dietary excess",
  "negative regulation of DNA binding",
  "response to zinc ion",
  "regulation of ATP biosynthetic process",
  "response to testosterone",
  "polyol catabolic process",
  "nucleoside triphosphate biosynthetic process",
  "specification of animal organ identity",
  "multi-multicellular organism process",
  "pattern specification process",
  "multi-organism reproductive process",
  "embryonic organ development",
  "response to ethanol",
  "stem cell development",
  "neural crest cell differentiation",
  "positive regulation of cholesterol efflux",
  "mitotic spindle assembly",
  "cerebral cortex development",
  "positive regulation of lipid localization",
  "actin filament capping",
  "negative regulation of protein depolymerization",
  "microtubule-based movement",
  "negative regulation of protein-containing complex disassembly",
  "protein polymerization",
  "regulation of cell shape",
  "signal transduction",
  "system process",
  "sensory perception",
  "glycosylation",
  "positive regulation of apoptotic process",
  "humoral immune response"
)

picked_tb3 <- tb3 %>% filter(description %in% functions_of_interest)
picked_tb4 <- tb4 %>% filter(description %in% functions_of_interest)
picked_tb5 <- tb5 %>% filter(description %in% functions_of_interest)
picked_tb6 <- tb6 %>% filter(description %in% functions_of_interest)
picked_tb7 <- tb7 %>% filter(description %in% functions_of_interest)
picked_tb8 <- tb8 %>% filter(description %in% functions_of_interest)
picked_tb9 <- tb9 %>% filter(description %in% functions_of_interest)

combined <- bind_rows(picked_tb9, picked_tb7, picked_tb6,  picked_tb8, picked_tb3, picked_tb5, picked_tb4 ) %>%
  filter(!is.na(p_adjust)) %>%
  mutate(
    source = factor(source, levels = c( "lesscode", "sameYR_antiWS", "antiYR_sameWS",  "antiYR_antiWS", "YRWS", "lowYR_sameWS", "sameYR_lowWS"  )),
    description = factor(description, levels = unique(description[order(source, p_adjust)]))
  )
ggplot(combined, aes(x = source, y = description)) +
  geom_point()
ggplot(combined, aes(x = source, y = description)) +
  geom_point(aes(size = fold_enrichment, fill = -log10(p_adjust)),
             shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient(low = "blue", high = "red", name = "-log10(p_adj)") +
  scale_size(range = c(1, 5)) +
  labs(x = NULL, y = NULL, title = "Nuc_type") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


### R

library(rGREAT)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(forcats)
setwd("/Path/to/Title/Library/E14")
list.files()


bed_WDR5_same <- read.table("./TSS_same_WDR5_M1.bed")
gr_WDR5_same <- GRanges(seqnames = bed_WDR5_same[, 1], ranges = IRanges(bed_WDR5_same[, 2], bed_WDR5_same[, 3]))
bed_bg <- read.table("./TSS_all_phase_adj+1Nuc_Di.bed")
gr_bg <- GRanges(seqnames = bed_bg[, 1], ranges = IRanges(bed_bg[, 2], bed_bg[, 3]))

res_WDR5_same <- great(gr_WDR5_same, "GO:BP", "Gencode_v38", background = gr_bg)
tb3 <- getEnrichmentTable(res_WDR5_same) %>% mutate(source = "WDR5_same")
write.csv(tb3, "tb3_results.csv", row.names = FALSE)

bed_WDR5_oppo <- read.table("./TSS_oppo_WDR5_M1.bed")
gr_WDR5_oppo <- GRanges(seqnames = bed_WDR5_oppo[, 1], ranges = IRanges(bed_WDR5_oppo[, 2], bed_WDR5_oppo[, 3]))

res_WDR5_oppo <- great(gr_WDR5_oppo, "GO:BP", "Gencode_v38", background = gr_bg)
tb4 <- getEnrichmentTable(res_WDR5_oppo) %>% mutate(source = "WDR5_oppo")
write.csv(tb4, "tb4_results.csv", row.names = FALSE)

tb4 <- read.csv("tb4_results.csv")
tb3 <- read.csv("tb3_results.csv")

functions_of_interest <- c(
  "cytoplasmic translation",
  "translation",
  "protein metabolic process",
  "ribonucleoprotein complex biogenesis",
  "ribosome biogenesis",
  "positive regulation of translation"      
)

picked_tb3 <- tb3 %>% filter(description %in% functions_of_interest)
picked_tb4 <- tb4 %>% filter(description %in% functions_of_interest)

combined <- bind_rows(picked_tb3, picked_tb4 ) %>%
  filter(!is.na(p_adjust)) %>%
  mutate(
    source = factor(source, levels = c("WDR5_same", "WDR5_oppo")),
    description = factor(description, levels = unique(description[order(source, p_adjust)]))
  )
ggplot(combined, aes(x = source, y = description)) +
  geom_point()
ggplot(combined, aes(x = source, y = description)) +
  geom_point(aes(size = fold_enrichment, fill = -log10(p_adjust)),
             shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient(low = "blue", high = "red", name = "-log10(p_adj)") +
  scale_size(range = c(1, 5)) +
  labs(x = NULL, y = NULL, title = "WDR5_associated_function") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

