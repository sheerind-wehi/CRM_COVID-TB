
setwd("/Users/sheerin.d/OneDrive - wehi.edu.au/2020/COVID-19/meta_analysis")

#load packages
library(gplots)
library(biomaRt)
library(edgeR)
library(org.Hs.eg.db)
library(pheatmap)
library(SummarizedExperiment)
library(TBSignatureProfiler)
library(curatedTBData)

#Wilk et. al.
#Single-cell
#PBMCs
#Severe - 4 ARDS, 4 less severely ill
#down in CD14 monocytes
wilk_cd14_mono_down <- c("TXNIP", "CD74", "HLA-DRB1", "HLA-DRA", "HLA-DPB1", "HLA-DRB5", "HLA-DMB", "EIF4B", "AP1S2", "IER2", "HLA-DQB1", "ZFP36L2", "TNFAIP2", "GNB2L1", "SLC25A6", "FOS", "EEF2", "PABPC1", "EEF1A1", "RGS2", "TPT1", "EIF1", "DUSP1", "JUNB")
#up in CD14 monocytes
wilk_cd14_mono_up <- c("IFI27", "S100A8", "MX2", "SIGLEC1", "XAF1", "ACSL1", "DYSF", "IFITM3", "MX1", "IFI6", "OAS3", "IFI44L", "OAS2", "PIM1", "VCAN", "RNASE2", "LGALS1", "PKM", "ACTB", "CLU", "PLAC8", "PLBD1", "SELL", "S100A12", "S100A9")
#down in NK cells
wilk_nk_down <- c("FCGR3A", "AHNAK", "FGFBP2", "GNLY", "JAK1", "GNB2L1", "ETS1", "SYNE1", "EEF2", "EIF4B", "ARL4C", "TPT1", "PABPC1", "C1orf63")
#up in NK cells
wilk_nk_up <- c("LAG3", "HAVCR2", "PLEK", "CD38", "HBB", "S100A9", "IGLC3", "GZMB", "ADAR", "MX2", "SP100", "MYOM2", "TRGC1", "IL32", "CD3D", "CD3G", "TRGC2", "RSAD2", "TRIM22", "PARP14", "CD38", "PARP9", "SAMD9", "SAMD9L", "EIF2AK2", "IFI6", "IFIT3", "OAS3", "OAS2", "IFI44", "ISG15", "CX3CR1", "ACTB", "AC009501.4", "XAF1", "IFI44L", "MX1", "RNA28S5", "STAT1")
#ISG signature
wilk_isg_sig <- c("IFI44", "IFI44L", "ISG15", "OAS2", "MX1", "IFIT3", "OAS3")
#up in activated granulocytes
wilk_act_gran <- c("ELANE", "LTF", "MMP8", "CEACAM8", "DEFA3", "LCN2")

#Huang et. al.
#Single-cell
#PBMCs
#1 critical, 1 severe, 5 moderate, 2 cured (plus other illnesses)
#IFN pathway genes
huang_ifn_path <- c("IFI27", "IFI44L", "IFI6", "IFIT2", "IFITM1", "IFITM2", "IFITM3", "IRF1", "IRF7", "ISG15", "ISG20", "MT2A", "OAS1", "SAMHD1", "ADAR", "BST2", "DDX50", "EPSTI1", "FCGR1A", "GBP1", "GBP2", "IFI35", "IFI44", "IFIT1", "IFIT3", "IRF2", "MX1", "MX2", "OAS2", "OAS3", "OASL", "PLSCR1", "RSAD2", "SP110", "TRIM22", "USP18", "XAF1")
#MAPK pathway genes
huang_mapk_path <- c("ATF3", "DDIT4", "DUSP1", "DUSP2", "FOS", "FOSB", "JUN", "JUNB", "JUND", "PRKCD")

#Wen et. al.
#Single-cell
#PBMCs
#10 recovering patients - early and late recovery stage
#up in monocytes
wen_mono <- c("IL1B", "JUN", "FOS", "JUNB", "KLF6", "CCL4", "CXCR4", "IRFD1", "IRF1", "IFI6", "IFITM3")
#up in NK cells
wen_nk <- c("NCAM1", "KLRF1", "KLRC1", "KLRD1", "CD56", "CD16")
#up in CD4 T cells
wen_cd4 <- c("CD3E", "CD4", "CCR7", "LEF1", "TCF7", "AQP3", "CD69", "CCR6", "CXCR6", "CCL5", "PRDM1", "FOS", "JUN", "KLF6", "S100A8", "IFITM3", "IFI6", "IL1B")
#up in CD8 T cells
wen_cd8 <- c("CCR7", "LEF1", "TCH7", "GZMK", "GZMB", "GNLY", "PRF1")
#up in B cells
wen_b <- c("CD19", "CD20", "IGHD", "IGHM", "IL4R", "TCL1A", "CD86", "CD27", "CD38", "IGHG", "XBP1", "MZB1")

#Liao et. al.
#Single-cell
#BALF
#3 severe patients, 3 mild
#group 1 macrophages (MDMs)
liao_g1_mac <- c("S100A9","S100A8","S100A12","VCAN","FCN1","CORO1A","SELL","CD14","CFP","RNASE2","SERPINB1","FPR1","COTL1","MPEG1","LST1","STAB1","RANSE6","MS4A6A")
#group 1/2 macrophages
liao_g1_2_mac <- c("IL1RN","CCL7","IFITM2","IFIT2","PLAC8","IFIT3","SERPINB9","NAMPT","LILRA5","IER2","ITITM3","MX1","NFKBIA","ANXA2R","RSAD2","IFIT1","DUSP6","ACTB")
#group 2 macrophages (pre-fibrotic)
liao_g2_mac <- c("SLC25A37","TNFSF10","CLU","CXCL10","APOBEC3A","CCL3L1","FFAR2","HSPA6","CCL2","HSPA1A","CCL8","DNAJB1","HSPA1B","CCL4","CCL3","CCL4L2","HSPB1","IDO1","ISG15","ISG20","TIMP1","NINJ1","VAMP5","CYP1B1","HSPH1","TYMP","CXCL11","WARS","BAG3","GBP1","SOD2","HSP90AA1","BCL2A1","GBP5","PLEK","NCF1","CALHM6","SLAMF7","SGK1","DEFB1","ANKRD22","GCH1","PLA2G7","CTSB")
#group 3 macrophages (fibrotic)
liao_g3_mac <- c("HAMP","LGMN","RGS1","SPP1","RNASE1","HMOX1","GPR183","ARL4C","C1QC","SDS","TGFB1","A2M","FPR3","PLTP","NR1H3","TTYH3","ABCA1","MARCKS","CTSZ","CD84","NRP2","CREG1","SDC3","CTSL","MMP14","SMPDL3A","CCL13","TMEM176B","TMEM176A","TREM2","LHFPL2","PMP22","PLD3","LIPA","MS4A4A","CSF1R","CD86","GPNMB","CREM","APOE","HSP90B1","C1QB","C1QC","CCL18","IL7R","CXCL8","IFITM1")
#group 4 macrophages (lung alveolar)
liao_g4_mac <- c("FABP5","FABP4","NURP1","APOC1","GCHFR","INHBA","ALOX5AP","ALDH2","CD52","HLA-DQA1","HLA-DQB1","LPL","MCEMP1","MARCO","TFRC","HPGD","FBP1","RBP4","ACP5","PEBP1","MSR1","MRC1","AKR1C3")
#down in expanded CD8 T cells in BALF
liao_exp_cd8_bal_down <- c("SELL", "ISG20", "CD27", "SAT1", "ISG15", "CD38", "CD74", "MT2A", "CD44")
#up in expanded CD8 T cells in BALF
liao_exp_cd8_bal_up <- c("GRP25", "KLRC1", "ZNF683", "CXCR6", "HOPX", "XCL1", "ITGAE", "ZYX")

#Xiong et. al.
#Bulk
#PBMCs
#3 patients, 3 healthy donors
#down in patients' PBMCs compared with healthy controls
xiong_pbmc_down <- c("GAA", "STT3A", "MGAT1", "LAIR1", "ADA2", "STEAP3", "PLA2G15", "CTSD", "PGD", "SLC1A5", "MARCO", "HK1", "C1QA", "MRC1", "C1QC")
#up in patients' PBMCs compared with healthy controls
xiong_pbmc_up <- c("MALAT1", "AC079793.1", "TSHZ2", "ZNF91", "SNORA73B", "RN7SL1", "FP671120.7", "RN7SL5P", "MT-ND6", "KIF21A", "IL5RA", "PRSS23", "CCDC80", "CCL4L2", "ZBTB20", "C1orf21", "SYTL2", "SLC4A4", "AGAP1", "AUTS2", "RHOBTB3", "CYBRD1", "TIMM8B", "CD180", "TNFSF10", "AK4", "MT1E", "IGLC2", "SDC1", "CKS2", "PRDM1", "TRIM7", "HTRA1", "GASK1B", "NT5DC2", "TMIGD3", "TPST1", "SRGAP1")
#down in patients' BALF compared with healthy controls
xiong_balf_down <- c("MALAT1", "AC079793.1", "TSHZ2", "ZNF91", "SNORA73B", "RN7SL1", "FP671120.7", "RN7SL5P", "MT-ND6", "KIF21A", "IL5RA", "PRSS23", "CCDC80", "CCL4L2", "ZBTB20", "C1orf21", "SYTL2", "SLC4A4", "AGAP1", "AUTS2", "RHOBTB3")
#up in patients' BALF compared with healthy controls
xiong_balf_up <- c("GAA", "STT3A", "MGAT1", "LAIR1", "ADA2", "STEAP3", "PLA2G15", "CTSD", "PGD", "SLC1A5", "MARCO", "HK1", "C1QA", "MRC1", "C1QC", "CYBRD1", "TIMM8B", "CD180", "TNFSF10", "AK4", "MT1E", "IGLC2", "SDC1", "CKS2", "PRDM1", "TRIM7", "HTRA1", "GASK1B", "NT5DC2", "TMIGD3", "TPST1", "SRGAP1")
#cytokines down in patients' PBMCs compared with healthy controls
xiong_cyto_pbmc_down <- c("TNFSF8", "FAM3B", "TNFSF13B", "CCL23", "CXCL3", "IL18", "TGFB1", "IL16", "GRN", "CXCL16", "BMP1", "GPI", "INHBA", "TNFSF12", "NRP1", "FAM3C", "CMTM3", "SECTM1", "IL1B", "IL1A", "HMGB1", "AIMP1", "CMTM7", "IL15")
#cytokines up in patients' PBMCs compared with healthy controls
xiong_cyto_pbmc_up <- c("CXCL11", "CXCL2", "TGFB2", "CCL3L1", "CXCL6", "CCL7", "CCL4", "CCL8", "SCGB3A1", "BMP3", "GDF15", "TNFSF10", "CCL2", "CXCL1", "IL33", "SPP1", "IL10", "IL36G", "TNFSF15", "IL36RN", "CCL5", "CXCL8", "CCL3", "NAMPT", "TNFSF14", "CCL18", "CXCL10", "VEGFA", "CCL20", "C5", "CXCL9", "TIMP1", "IL1RN")
#cytokines down in patients' BALF compared with healthy controls
xiong_cyto_balf_down <- c("GDF7", "CCL28", "IL24", "CD40LG", "GDF11", "LTB", "IL23A", "FLT3LG", "LTA", "IL32", "FASLG", "WNT7A", "NOG", "XCL2", "XCL1", "CSF1", "IL1B", "OSM", "CCL4", "NAMPT", "CXCL1", "CMTM2", "CXCL8", "IL1RN", "CCL3L1", "TGFB3", "CLCF1", "IFNG", "WNT1", "CCL5", "IL7", "CXCL5")
#cytokines up in patients' BALF compared with healthy controls
xiong_cyto_balf_up <- c("IL10", "CXCL10", "CXCL9", "TNFSF10", "IL15", "NRG1", "TIMP1", "TNFSF13", "GRN", "GPI", "BMP8A", "TNFSF8", "CMTM7", "PF4", "C5", "IL18", "EBI3", "CMTM5", "AREG", "VSTM1", "CXCL2", "CXCL16", "BMP6", "TNFSF14", "CCL3", "IL6", "NRP1", "VEGFA", "SECTM1", "BMP1", "NODAL", "CMTM3", "IL12A", "TNFSF13B", "IL16", "THNSL2", "TNF", "CNTF", "PF4V1", "PPBP", "TNFSF4", "TGFB1", "TNFSF12", "BMP8B", "FAM3C", "HMGB1", "AIMP1")

#Hadjadj et. al.
#Nanostring
#Whole blood
#18 critical, 17 severe, 15 mild to moderate
#down in patients' whole blood
hadjadj_wb_down <- c("MX1", "IFITM1", "IFIT2")
#up in patients' whole blood
hadjadj_wb_up <- c("IFNAR1", "JAK1", "TYK2", "IL6R", "SOCS3", "STAT3", "TNFSF10", "IL1B", "IL1RA", "IL10", "CXCR2", "IRF1", "PSMB8", "IFNAR2", "IFNA2", "IRF4", "IRF5", "SOCS1", "STAT2", "IFI35", "IRF7", "BST2", "STAT1")
#interferon-stimulated gene signature
hadjadj_imp_isgs <- c("IFI44L", "IFI27", "RSAD2", "SIGLEC1", "IFIT1", "IS15")
#down in mild to moderate patients compared with healthy controls
hadjadj_mild_mod_down <- c("IL1RAP", "KLRB1", "SKI", "FCER1A")
#up in mild to moderate patients compared with healthy controls
hadjadj_mild_mod_up <- c("BST2", "CD36", "CASP1", "BST1", "JAK2", "TLR7", "IRF4", "TNFSF15", "CD14", "CDKN1A", "SOCS1", "FCGR1A", "IL1RN", "CCL2", "S100A8", "LAMP3", "TLR5", "PML", "GBP1", "CXCL10", "STAT2", "IFIT2", "SERPING1", "MX1", "C2", "IF135", "C1QB", "S100A9", "LILRA6", "TNFS10", "CCR1", "CEACAM1", "CLEC5A", "CDKN1A", "CD274")
#down in severe patients compared with mild to moderate patients
hadjadj_mod_v_sev_down <- c("CCR5", "MSR1", "TGFBI", "LAG3", "TLR7")
#up in severe patients compared with mild to moderate patients
hadjadj_mod_v_sev_up <- c("BTK", "CD46", "TLR9", "SELPLG", "TLR8", "CD82", "IL4R", "IL1R1", "TLR5", "IL18RAP", "CR1", "LY96", "ARG1", "S100A8", "IL32", "NCF4", "CD45R0", "LILRA5", "S100A9", "MME", "IFNAR1", "IL1R2", "IRAK3", "PPBP", "MAPK14")
#down in critical patients compared with severe patients
hadjadj_sev_v_crit_down <- c("GZMB", "KIR_ACTIVATING_SUBGROUP_1", "PRF1")
#up in critical patients compared with severe patients
hadjadj_sev_v_crit_up <- c("ITGB1", "CD55", "NOD2", "LGALS3", "MRC1", "ITLN1", "GP1BB", "CTSG", "AHR", "IL1RA", "CD9", "TAL1", "PPBP", "ITGA2B", "IL1RAP")

#Wei et. al.
#Single-cell
#PBMCs
#4 patients before, during, and after ICU care
#upregulated in inactivated monocytes
wei_inact_mono_up <- c("LYZ", "CST3", "S100A9", "S100A8", "IL8", "NEAT1", "HBA1", "HBA2", "HBB", "IFIT3", "IFITM1", "IRF7", "OAS1", "OAS2", "OAS3", "RSAD2")
#upregulated in classical monocytes
wei_class_mono_up <- c("RP11-1143G9.4", "LYZ", "CST3", "S100A9", "S100A8", "IL8", "NEAT1", "CMPK2", "COPB2", "DDX58", "IFIT1", "IFIT3", "IFITM1", "IRF7", "OAS1", "OAS2", "OAS3", "RSAD2", "USP18")
#upregulated in B cells
wei_b_up <- c("MS4A1", "BANK1", "RALGPS2", "CD79A", "LINC00926", "HLA-DQA1", "MZB1", "IGJ", "IGLL5")
#upregulated in T cells
wei_t_up <- c("IL7R", "RORA", "IL32", "LEF1", "BCL11B", "GNLY", "CCL5", "NKG7")


#Lee et. al.
#Single-cell
#PBMCs
#8 COVID-19 patients (severe, mild, asymptomatic), 5 influenza, 4 HCs
#COVID-19 vs. HC
lee_cov_v_hc <- read.csv("COVID_v_HC.csv")
lee_cov_v_hc <- as.character(lee_cov_v_hc$Gene)
#COVID-19-specific
COVID_spec <- read.csv("COVID_specific.csv")
#upregulated in COVID-19 only
lee_cov_spec_up <- COVID_spec[-grep("-1", COVID_spec$Direction),]
lee_cov_spec_up <- as.character(lee_cov_spec_up$Gene)
#downregulated in COVID-19 only
lee_cov_spec_down <- COVID_spec[grep("-1", COVID_spec$Direction),]
lee_cov_spec_down <- as.character(lee_cov_spec_down$Gene)
#flu vs. HC
#COVID-19 vs. HC
lee_flu_v_hc <- read.csv("Flu_v_HC.csv")
lee_flu_v_hc <- as.character(lee_flu_v_hc$Gene)
#flu-specific
FLU_spec <- read.csv("Flu_specific.csv")
#upregulated in flu only
lee_flu_spec_up <- FLU_spec[-grep("-1", FLU_spec$Direction),]
lee_flu_spec_up <- as.character(lee_flu_spec_up$Gene)
#downregulated in COVID-19 only
lee_flu_spec_down <- FLU_spec[grep("-1", FLU_spec$Direction),]
lee_flu_spec_down <- as.character(lee_flu_spec_down$Gene)

#Lee et. al.
#Single-cell
#PBMCs
#blood samples from 3 control patients with flu-like symptoms (SARS-CoV-2 negative), and 3 SARS2-CoV-2  positive  patients
#1  outpatient  with  mild  disease  and  2 patients  with  severe disease  admitted  to  ICU
#whole blood
silvin_comb <- read.csv("Silvin_comb.csv", row.names = 1, header = T)
#up-regulated
silvin_comb_up <- silvin_comb[silvin_comb$logFC > 0,]
var_silvin_comb <- apply(silvin_comb_up, 1, var)
select_var_silvin_comb <- names(sort(var_silvin_comb, decreasing = T))[1:200]
silvin_comb_up <- silvin_comb_up[select_var_silvin_comb,]
silvin_comb_up <- as.character(rownames(silvin_comb_up))
#down-regulated
silvin_comb_down <- silvin_comb[silvin_comb$logFC < 0,]
silvin_comb_down <- silvin_comb_down[,1]
#monocytes
silvin_mono <- read.csv("Silvin_mono.csv", header = T)
#up-regulated
silvin_mono_up <- silvin_mono[silvin_mono$avg_logFC > 0,]
silvin_mono_up <- silvin_mono_up[,1]
#down-regulated
silvin_mono_down <- silvin_mono[silvin_mono$avg_logFC < 0,]
silvin_mono_down <- silvin_mono_down[,1]
#neutrophils
silvin_neut <- read.csv("Silvin_neut.csv", header = T)
#up-regulated
silvin_neut_up <- silvin_neut[silvin_neut$avg_logFC > 0,]
silvin_neut_up <- silvin_neut_up[,1]
#down-regulated
silvin_neut_down <- silvin_neut[silvin_neut$avg_logFC < 0,]
silvin_neut_down <- silvin_neut_down[,1]

#Pulendran et. al.
#RNA-seq
#Whole blood
bali_btms <- c("TAP1","IRF7","USP18","HESX1","OAS1","EIF2AK2","HERC5","IFIT2","SERPING1","BCL2A1","IFNGR2","MGST1","APOBEC3A","STAT1","IFI27","DDX60","SP100","RSAD2","IFIT3","IFIT1","C1QB","OASL","SIGLEC1","OAS3","HLX","DDX58","IFIH1","MX2","PML","PLSCR1","PARP9","IFITM1","ATF3")
bali_cov <- read.csv("sig_res_COVID_v_healthy_BALI_bulk_PBMCs.csv", row.names = 1)
var_bali_cov <- apply(bali_cov, 1, var)
select_var_bali_cov <- names(sort(var_bali_cov, decreasing = T))[1:200]
bali_cov <- bali_cov[select_var_bali_cov,]
bali_cov <- as.character(rownames(bali_cov))
bali_mod <- read.csv("sig_res_MOD_BALI_bulk_PBMCs.csv", row.names = 1)
var_bali_mod <- apply(bali_mod, 1, var)
select_var_bali_mod <- names(sort(var_bali_mod, decreasing = T))[1:200]
bali_mod <- bali_mod[select_var_bali_mod,]
bali_mod <- as.character(rownames(bali_mod))
bali_sev <- read.csv("sig_res_SEV_BALI_bulk_PBMCs.csv", row.names = 1)
var_bali_sev <- apply(bali_sev, 1, var)
select_var_bali_sev <- names(sort(var_bali_sev, decreasing = T))[1:200]
bali_sev <- bali_sev[select_var_bali_sev,]
bali_sev <- as.character(rownames(bali_sev))
bali_icu <- read.csv("sig_res_ICU_BALI_bulk_PBMCs.csv", row.names = 1)
var_bali_icu <- apply(bali_icu, 1, var)
select_var_bali_icu <- names(sort(var_bali_icu, decreasing = T))[1:200]
bali_icu <- bali_icu[select_var_bali_icu,]
bali_icu <- as.character(rownames(bali_icu))

#Guo et al.
guo_mono <- read.csv("guo_mono.csv")
guo_sev_mono <- as.character(guo_mono[2:nrow(guo_mono),1])[1:40]
guo_rem_hc_mono <- as.character(guo_mono[2:789,4])[1:100]
guo_T <- read.csv("guo_t.csv")
guo_sev_T <- as.character(guo_T[2:583,1])[1:40]
guo_rem_hc_T <- as.character(guo_T[2:374,3])[1:100]

#O'Garra flu data
og_flu <- read.csv("ogarra_flu_WB.csv", row.names = 1, header = T)
var_og_flu <- apply(og_flu, 1, var)
select_var_og_flu <- names(sort(var_og_flu, decreasing = T))[1:200]
og_flu <- og_flu[select_var_og_flu,]
og_flu <- as.character(rownames(og_flu))

#Kaforou 3-gene
gene3 <- c("HERC6", "IGF1R", "NAGK")


#load O'Garra et. al. TB data
ogarra_SE <- readRDS(file = "/Users/sheerin.d/OneDrive - wehi.edu.au/2020/Bioinformatics/HHC_study/kristina_analysis/ogarra_SE.rds")
OG_genes <- assays(ogarra_SE)[[4]]
keep <- rowSums(OG_genes > 2) >= 5
summary(keep)
#   Mode   FALSE    TRUE 
#logical  12163   13206 
OG_genes <- OG_genes[keep,]
OG_var <- apply(OG_genes, 1, var)
select_var_5000 <- names(sort(OG_var, decreasing = T))[1:5000]
OG_genes <- OG_genes[select_var_5000,]
write.csv(rownames(OG_genes), "OG_genes_5000.csv")


ogarra_SE$TBStatus[ogarra_SE$TBStatus == "Leicester_Control"] = "3. Leicester Contact"
ogarra_SE$TBStatus[ogarra_SE$TBStatus == "Leicester_Latent"] = "4. Leicester Latent"
ogarra_SE$TBStatus[ogarra_SE$TBStatus == "Leicester_LTBI"] = "4. Leicester Latent"
ogarra_SE$TBStatus[ogarra_SE$TBStatus == "London_Latent"] = "1. London Latent"
ogarra_SE$TBStatus[ogarra_SE$TBStatus == "SA_Latent"] = "2. South Africa Latent"
ogarra_SE$TBStatus[ogarra_SE$TBStatus == "Leicester_Progressor"] = "5. Leicester Progressor"
ogarra_SE$TBStatus[ogarra_SE$TBStatus == "London_PTB"] = "6. London PTB"
ogarra_SE$TBStatus[ogarra_SE$TBStatus == "SA_PTB"] = "7. South Africa PTB"
ogarra_SE$TBStatus[ogarra_SE$TBStatus == "Active_TB"] = "8. Active TB"
TBStatus <- ogarra_SE$TBStatus

#compile COVID-19 signatures
COVID_sigs <- list(Wilk_sc_PBMC_monocytes_up = wilk_cd14_mono_up, Wilk_sc_PBMC_monocytes_down = wilk_cd14_mono_down, Wilk_sc_PBMC_NK_cells_up = wilk_nk_up, Wilk_sc_PBMC_NK_cells_down = wilk_nk_down, Wilk_sc_PBMCs_ISG_signature = wilk_isg_sig, Wilk_sc_PBMC_activated_granulocytes = wilk_act_gran, Huang_sc_PBMC_IFN_signature = huang_ifn_path, Wen_sc_PBMC_monocytes = wen_mono, Wen_sc_PBMC_NK_cells = wen_nk, Wen_sc_PBMC_CD4_T_cells = wen_cd4, Wen_sc_PBMC_CD8_T_cells = wen_cd8, Wen_sc_PBMC_B_cells = wen_b, Xiong_bulk_PBMC_gene_signature_up = xiong_balf_up, Xiong_bulk_PBMC_gene_signature_down = xiong_balf_down, Xiong_sc_PBMC_cytokines_up = xiong_cyto_balf_up, Xiong_sc_PBMC_cytokines_down = xiong_cyto_balf_down, Liao_sc_BALF_G1_macrophages = liao_g1_mac, Liao_sc_BALF_G1_2_macrophages = liao_g1_2_mac, Liao_sc_BALF_G2_macrophages = liao_g2_mac, Liao_sc_BALF_G3_macrophages = liao_g3_mac, Liao_sc_BALF_G4_macrophages = liao_g4_mac, Liao_sc_BALF_CD8_T_cells_up = liao_exp_cd8_bal_up, Liao_sc_BALF_CD8_T_cells_down = liao_exp_cd8_bal_down, Hadjadj_nanostring_WB_gene_signature_up = hadjadj_wb_up, Hadjadj_nanostring_WB_gene_signature_down = hadjadj_wb_down, Hadjadj_nanostring_WB_ISG_signature = hadjadj_imp_isgs, Hadjadj_nanostring_WB_mild_moderate_up = hadjadj_mild_mod_up, Hadjadj_nanostring_WB_mild_moderate_down = hadjadj_mild_mod_down, Hadjadj_nanostring_WB_severe_up = hadjadj_mod_v_sev_up, Hadjadj_nanostring_WB_severe_down = hadjadj_mod_v_sev_down, Hadjadj_nanostring_critical_up = hadjadj_sev_v_crit_up, Hadjadj_nanostring_critical_down = hadjadj_sev_v_crit_down, Wei_sc_PBMC_inactivated_monocytes = wei_inact_mono_up, Wei_sc_PBMC_classical_monocytes = wei_class_mono_up, Wei_sc_PBMCs_T_cells = wei_t_up, Wei_sc_PBMC_B_cells = wei_b_up, Silvin_sc_WB_combined_signature = silvin_comb_up, Silvin_sc_WB_monocytes_up = silvin_mono_up, Silvin_sc_WB_monocytes_down = silvin_mono_down, Silvin_sc_WB_neutrophils_up = silvin_neut_up, Silvin_sc_WB_neutrophils_down = silvin_neut_down, Arunachalam_bulk_PBMC_blood_modules = bali_btms, Arunachalam_bulk_PBMC_covid_combined = bali_cov, Arunachalam_bulk_PBMC_moderate = bali_mod, Arunachalam_bulk_PBMC_severe = bali_sev, Arunachalam_bulk_PBMC_intensive_care = bali_icu, Guo_sc_PBMC_monocytes_severe = guo_sev_mono, Guo_sc_PBMC_monocytes_healthy = guo_rem_hc_mono, Guo_sc_PBMC_T_cells_severe = guo_sev_T, Guo_sc_PBMC_T_cells_healthy = guo_rem_hc_T, Dunning_bulk_WB_flu = og_flu)

COVID_rearr <- list(Wilk_et_al._PBMCs_monocytes = wilk_cd14_mono_up, Wen_et_al._PBMCs_monocytes = wen_mono, Wei_et_al._PBMCs_inactivated_monocytes = wei_inact_mono_up, Wei_et_al._PBMCs_classical_monocytes = wei_class_mono_up, Silvin_et_al._WB_monocytes = silvin_mono_up, Wilk_et_al._PBMCs_NK_cells = wilk_nk_up, Liao_et_al._BALF_G1_macrophages = liao_g1_mac, Liao_et_al._BALF_G1_2_macrophages = liao_g1_2_mac, Liao_et_al._BALF_G2_macrophages = liao_g2_mac, Liao_et_al._BALF_G3_macrophages = liao_g3_mac, Liao_et_al._BALF_G4_macrophages = liao_g4_mac, Wen_et_al._PBMCs_NK_cells = wen_nk, Wen_et_al._PBMCs_CD4_T_cells = wen_cd4, Wen_et_al._PBMCs_CD8_T_cells = wen_cd8, Liao_et_al._BALF_CD8_T_cells = liao_exp_cd8_bal_up, Wen_et_al._T_cells = guo_sev_T, Wen_et_al._PBMCs_B_cells = wen_b, Wei_et_al._PBMCs_B_cells = wei_b_up, Wilk_et_al._PBMCs_activated_granulocytes = wilk_act_gran, Silvin_et_al._WB_neutrophils = silvin_neut_up, Silvin_et_al._WB_combined_signature = silvin_comb_up, Wilk_et_al._PBMCs_ISG_signature = wilk_isg_sig, Huang_et_al._PBMCs_IFN_signature = huang_ifn_path, Hadjadj_et_al._WB_ISG_signature = hadjadj_imp_isgs, Hadjadj_et_al._WB_gene_signature = hadjadj_wb_up, Hadjadj_et_al._WB_mild_moderate = hadjadj_mild_mod_up, Hadjadj_et_al._WB_severe = hadjadj_mod_v_sev_up, Hadjadj_et_al._WB_critical = hadjadj_sev_v_crit_up, Xiong_et_al._PBMCs_gene_signature = xiong_balf_up, Xiong_et_al._PBMCs_cytokines = xiong_cyto_balf_up, bali_cov = bali_cov, bali_mod = bali_mod, bali_sev = bali_sev, bali_icu = bali_icu, bali_btms = bali_btms, og_flu = og_flu)

saveRDS(COVID_sigs, "COVIDsignatures.rds")

library(data.table)
lapply(COVID_rearr, function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))

#run TBSignatureProfiler
og_scores_GSVA <- runTBsigProfiler(ogarra_SE, useAssay = "log_cpm", signatures = COVID_rearr, algorithm = "GSVA")
og_scores_ssGSEA = runTBsigProfiler(ogarra_SE, useAssay = "log_cpm", signatures = COVID_sigs, algorithm = "ssGSEA")

SA_latent_pval <- tableAUC(og_scores_GSVA[,grep("1. London Latent|2. South Africa Latent", og_scores_GSVA$TBStatus)],
                             annotationColName = "TBStatus",
                             signatureColNames = names(COVID_rearr),
                             num.boot = 100,
                             pb.show = FALSE,
                             output = "data.frame")
write.csv(SA_latent_pval, "SA_latent_pval.csv")
leic_contact_pval <- tableAUC(og_scores_GSVA[,grep("1. London Latent|3. Leicester Contact", og_scores_GSVA$TBStatus)],
                             annotationColName = "TBStatus",
                             signatureColNames = names(COVID_rearr),
                             num.boot = 100,
                             pb.show = FALSE,
                             output = "data.frame")
write.csv(leic_contact_pval, "leic_contact_pval.csv")
leic_latent_pval <- tableAUC(og_scores_GSVA[,grep("1. London Latent|4. Leicester Latent", og_scores_GSVA$TBStatus)],
                            annotationColName = "TBStatus",
                            signatureColNames = names(COVID_rearr),
                            num.boot = 100,
                            pb.show = FALSE,
                            output = "data.frame")
write.csv(leic_latent_pval, "leic_latent_pval.csv")
leic_prog_pval <- tableAUC(og_scores_GSVA[,grep("1. London Latent|5. Leicester Progressor", og_scores_GSVA$TBStatus)],
                           annotationColName = "TBStatus",
                           signatureColNames = names(COVID_rearr),
                           num.boot = 100,
                           pb.show = FALSE,
                           output = "data.frame")
write.csv(leic_prog_pval, "leic_prog_pval.csv")
lon_PTB_pval <- tableAUC(og_scores_GSVA[,grep("1. London Latent|6. London PTB", og_scores_GSVA$TBStatus)],
                           annotationColName = "TBStatus",
                           signatureColNames = names(COVID_rearr),
                           num.boot = 100,
                           pb.show = FALSE,
                           output = "data.frame")
write.csv(lon_PTB_pval, "lon_PTB_pval.csv")
SA_PTB_pval <- tableAUC(og_scores_GSVA[,grep("1. London Latent|7. South Africa PTB", og_scores_GSVA$TBStatus)],
                         annotationColName = "TBStatus",
                         signatureColNames = names(COVID_rearr),
                         num.boot = 100,
                         pb.show = FALSE,
                         output = "data.frame")
write.csv(SA_PTB_pval, "SA_PTB_pval.csv")
leic_act_pval <- tableAUC(og_scores_GSVA[,grep("1. London Latent|8. Active TB", og_scores_GSVA$TBStatus)],
                        annotationColName = "TBStatus",
                        signatureColNames = names(COVID_rearr),
                        num.boot = 100,
                        pb.show = FALSE,
                        output = "data.frame")
write.csv(leic_act_pval, "leic_act_pval.csv")

og_GSVA_res <- as.data.frame(colData(og_scores_GSVA)[,c(4,23:58)])
og_GSVA_res <- reshape2::melt(og_GSVA_res, id.vars = c("TBStatus"))
levels(og_GSVA_res$variable) <- c("Wilk Monocytes", "Wen Monocytes", "Wen Inactivated Monocytes", "Wei Classical Monocytes", "Silvin Monocytes (WB)", "Wilk NK Cells", "Liao G1 M  (BALF)", "Liao G1/2 M  (BALF)", "Liao G2 M  (BALF)", "Liao G3 M  (BALF)", "Liao G4 M  (BALF)", "Wen NK cells", "Wen CD4 T cells", "Wen CD8 T cells", "Liao CD8 T cells", "Wen T cells", "Wen B Cells", "Wei B Cells", "Wilk Activated Granulocytes", "Silvin Neutrophils (WB)", "Silvin WB Signature", "Wilk ISG Signature", "Huang IFN Signature", "Hadjadj ISG Signature (WB)", "Hadjadj WB Signature", "Hadjadj Mild/Moderate (WB)", "Hadjadj Severe (WB)", "Hadjadj Critical (WB)", "Xiong PBMC Signature", "Xiong Cytokine Signature", "Arunachalam PBMC Signature", "Arunachalam Moderate", "Arunachalam Severe", "Arunachalam ICU", "Arunachalam BTMs", "Dunning Influenza (WB)")
ggplot(og_GSVA_res, aes(x = TBStatus, y = value, group = TBStatus)) +
  geom_boxplot(aes(fill = TBStatus), outlier.shape = NA) +
  geom_jitter(aes(), width = 0.2, size = 0.5) +
  facet_wrap(~variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  scale_fill_discrete(name = "Group") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "1. London Latent", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

library(rstatix)
stat.test <- og_GSVA_res %>%
  group_by(variable) %>%
  t_test(value ~ TBStatus, ref.group = "1. London Latent") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- as.data.frame(stat.test)
stat.test
write.csv(stat.test, "o_garra_COVID_stats2.csv")

signatureBoxplot(inputData = og_scores_GSVA, 
                 name = "Gene Set Variation Analysis of O'Garra et. al. TB data",
                 signatureColNames = names(COVID_rearr),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

signatureBoxplot(inputData = og_scores_ssGSEA, 
                 name = "ssGSEA",
                 signatureColNames = names(COVID_sigs),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#check 3-gene signature
gene3 <- list(gene3 = gene3)
gene3_scores_GSVA <- runTBsigProfiler(ogarra_SE, useAssay = "log_cpm", signatures = gene3, algorithm = "GSVA")
signatureBoxplot(inputData = gene3_scores_GSVA, 
                 name = "Gene Set Variation Analysis of Kaforou 3-gene COVID signature on TB data",
                 signatureColNames = names(gene3),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#assess on other datasets
library(curatedTBData)
library(MultiAssayExperiment)
library(BiocParallel)
library(dplyr)

data("DataSummary",package="curatedTBData")

geo_rnaseq <- DataSummary %>% dplyr::as_tibble() %>%
  dplyr::filter(stringr::str_detect(.data$GeneralType, "Illumina RNA-seq")) %>%
  dplyr::select(.data$`GEOAccession`)

obj_list_rseq <- curatedTBData(geo_access = geo_rnaseq$GEOAccession[c(1,2,7)]) 

param <- SerialParam(progressbar=TRUE)

obj_norm <- bplapply(obj_list_rseq, function(x)
  Normalization(x, microarray_method = "quantile", RNAseq_method = "TMM",
                experiment_name = "assay_raw"), BPPARAM = param) # Takes some time

obj_norm <- bplapply(obj_norm, function(x)
  x[,colData(x)[,"TBStatus"]!= "NA"], BPPARAM = param)

object_match <- bplapply(obj_norm, function(x)
  MatchProbe(x, UseAssay = c("TMM","quantile","RMA"),
             createExperimentName = "assay_MatchProbe"), BPPARAM = param)

multi_set_PTB_Latent <- bplapply(object_match, function(x)
  subset_curatedTBData(x, annotationColName = "TBStatus",
                       annotationCondition = c("Latent","PTB"),
                       experiment_name = "assay_MatchProbe"),
  BPPARAM = param) %>% plyr::compact()

GC6_meta <- as.data.frame(colData(multi_set_PTB_Latent$GSE94438))
write.csv(GC6_meta, "GC6_meta.csv")
ACS_meta <- as.data.frame(colData(multi_set_PTB_Latent$GSE79362))
write.csv(ACS_meta, "ACS_meta.csv")
OG_meta <- as.data.frame(colData(multi_set_PTB_Latent$GSE107994))
write.csv(OG_meta, "OG_meta.csv")

gsva_PTB_Latent <- lapply(multi_set_PTB_Latent,
                            function(x) TBSignatureProfiler::runTBsigProfiler(
                              input = x,
                              useAssay = assayNames(x),
                              signatures = COVID_rearr,
                              algorithm = "GSVA",
                              combineSigAndAlgorithm = TRUE))

library(ggplot2)
library(ggpubr)
#look at all the Leicester groups first
colnames(colData(gsva_PTB_Latent$GSE107994))
OG_GSVA <- as.data.frame(colData(gsva_PTB_Latent$GSE107994)[,c(4,22,26:62)])
OG_GSVA$Guo_et_al._PBMC_monocytes <- -(OG_GSVA$Guo_et_al._PBMC_monocytes)
levels(OG_GSVA$Timepoint_month) <- c("Non-progressor","Progressor <0.5","Non-progressor","Progressor 0.6-2.5","Non-progressor","Progressor 2.6-3.3","Progressor 3.4-5","Baseline")
levels(OG_GSVA$TBType) <- c("Contact (non-pulmonary)","Contact (pulmonary)","TB (non-pulmonary)","TB (pulmonary)")
OG_GSVA <- reshape2::melt(OG_GSVA, id.vars = c("TBStatus","TBType","Timepoint_month"))
levels(OG_GSVA$variable) <- c("Wilk Monocytes", "Wen Monocytes", "Wen Inactivated Monocytes", "Wei Classical Monocytes", "Silvin Monocytes (WB)", "Guo Monocytes", "Liao G1 M  (BALF)", "Liao G2 M  (BALF)", "Liao G3 M  (BALF)", "Liao G4 M  (BALF)", "Wilk NK Cells", "Wen NK cells", "Wen CD4 T cells", "Wen CD8 T cells", "Liao CD8 T cells", "Guo CD8 T cells", "Wen B Cells", "Wei B Cells", "Wilk Activated Granulocytes", "Silvin Neutrophils (WB)", "Silvin WB Signature", "Wilk ISG Signature", "Huang IFN Signature", "Hadjadj ISG Signature (WB)", "Hadjadj WB Signature", "Hadjadj Mild/Moderate (WB)", "Hadjadj Severe (WB)", "Hadjadj Critical (WB)", "Xiong PBMC Signature", "Xiong Cytokine Signature", "Arunachalam PBMC Signature", "Arunachalam Moderate", "Arunachalam Severe", "Arunachalam ICU", "Arunachalam BTMs", "Dunning Influenza (WB)")
ggplot(OG_GSVA, aes(x = TBType, y = value, group = TBType)) +
  geom_boxplot(aes(fill = TBType), outlier.shape = NA) +
  geom_jitter(aes(colour = Timepoint_month), width = 0.2) +
  scale_colour_manual(values = c("black","red","orange","darkgoldenrod1","brown","grey")) +
  facet_wrap(~ variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Contact (non-pulmonary)", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

stat.test01 <- OG_GSVA %>%
  group_by(variable) %>%
  t_test(value ~ TBType, ref.group = "Contact (non-pulmonary)") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test01 <- as.data.frame(stat.test01)
write.csv(stat.test01, "supp1_COVID_stats.csv")

colnames(colData(gsva_PTB_Latent$GSE94438))
GC6_GSVA <- as.data.frame(colData(gsva_PTB_Latent$GSE94438)[,c(4,22,24:29,43,59,45:50,58)])
GC6_GSVA <- reshape2::melt(GC6_GSVA, id.vars = c("TBStatus","TimeToTB_Month"))
GC6_GSVA$TimeToTB_Month <- as.numeric(GC6_GSVA$TimeToTB_Month)
GC6_GSVA$TimeToTB_Month[is.na(GC6_GSVA$TimeToTB_Month)] <- 20
levels(GC6_GSVA$variable) <- c("Wilk Monocytes","Wen Monocytes","Wei Inactivated Monocytes","Wei Classical Monocytes","Silvin Monocytes (WB)","Guo Monocytes","Silvin Neutrophils (WB)","Dunning Influenza (WB)","Wilk ISG Signature","Huang IFN Signature","Hadjadj ISG Signature","Hadjadj WB Signature","Hadjadj Mild/Moderate (WB)","Hadjadj Severe (WB)","Arunachalam BTMs")
ggplot(GC6_GSVA, aes(x = TBStatus, y = value, group = TBStatus)) +
  geom_boxplot(aes(fill = TBStatus), outlier.shape = NA) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  geom_jitter(aes(colour = TimeToTB_Month, size = -TimeToTB_Month), width = 0.2) +
  scale_colour_gradient2(low = "red", mid = "blue",
                         high = "grey", space = "Lab", midpoint = 10, na.value = "grey50", guide = "colourbar") +
  scale_size_continuous(range = c(0.5,3),) + 
  facet_wrap(~ variable, ncol = 7) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Latent", hide.ns = T, label.y = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

stat.test2 <- GC6_GSVA %>%
  group_by(variable) %>%
  t_test(value ~ TBStatus, ref.group = "Latent") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test2 <- as.data.frame(stat.test2)
write.csv(stat.test2, "GC6_COVID_stats.csv")

colnames(colData(gsva_PTB_Latent$GSE79362))
ACS_GSVA <- as.data.frame(colData(gsva_PTB_Latent$GSE79362)[,c(4,25,28:32,47,63,49:54,62)])
ACS_GSVA <- reshape2::melt(ACS_GSVA, id.vars = c("TBStatus","TimeToTB_Day"))
ACS_GSVA$TimeToTB_Day <- log2(ACS_GSVA$TimeToTB_Day)
PTB <- ACS_GSVA[grep("PTB", ACS_GSVA$TBStatus),]
PTB$TimeToTB_Day[is.na(PTB$TimeToTB_Day)] <- 10
ACS_GSVA[grep("PTB", ACS_GSVA$TBStatus),] <- PTB
ACS_GSVA$TimeToTB_Day[ACS_GSVA$TimeToTB_Day == "NaN"] <- 10
ACS_GSVA$TimeToTB_Day[is.na(ACS_GSVA$TimeToTB_Day)] <- 15
ACS_GSVA$TimeToTB_Day[ACS_GSVA$TimeToTB_Day == "-Inf"] <- 0
levels(ACS_GSVA$variable) <- c("Wilk Monocytes","Wen Monocytes","Wei Inactivated Monocytes","Wei Classical Monocytes","Silvin Monocytes (WB)","Silvin Neutrophils (WB)","Dunning Influenza (WB)","Wilk ISG Signature","Huang IFN Signature","Hadjadj ISG Signature","Hadjadj WB Signature","Hadjadj Mild/Moderate (WB)","Hadjadj Severe (WB)","Arunachalam BTMs")
ggplot(ACS_GSVA, aes(x = TBStatus, y = value, group = TBStatus)) +
  geom_boxplot(aes(fill = TBStatus), outlier.shape = NA) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  geom_jitter(aes(colour = TimeToTB_Day, size = -TimeToTB_Day), width = 0.2) +
  scale_colour_gradient2(low = "red", mid = "blue",
                         high = "grey", space = "Lab", midpoint = 7, na.value = "grey50", guide = "colourbar") +
  scale_size_continuous(range = c(0.5,3),) + 
  facet_wrap(~ variable, ncol = 7) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Latent", hide.ns = T, label.y = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

stat.test3 <- ACS_GSVA %>%
  group_by(variable) %>%
  t_test(value ~ TBStatus, ref.group = "Latent") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test3 <- as.data.frame(stat.test3)
write.csv(stat.test3, "ACS_COVID_stats.csv")

signatureBoxplot(inputData = gsva_PTB_Latent$GSE79362, 
                 name = "Gene Set Variation Analysis of Adolescent Cohort Study TB Data",
                 signatureColNames = names(COVID_sigs),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

signatureBoxplot(inputData = gsva_PTB_Latent$GSE107994, 
                 name = "Gene Set Variation Analysis of Grand Challenges 6 TB Data",
                 signatureColNames = names(COVID_sigs),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

get_auc_distribution(ssgsea_PTB_Latent_combine) + 
  ggtitle("Ridge plot of AUC for PTB vs. Latent") + 
  theme(axis.text.x = element_text(colour="Black", size=12, hjust = 0.5, 
                                   face="bold"),
        axis.text.y = element_text(size=12, angle = 0, hjust = 0.5))

#get O'Garra Leicester data
leic_TB <- get_curatedTBData(geo_access = "GSE107994", include.SCAN = FALSE)
leic_TB <- leic_TB$GSE107994
leic_raw <- assays(leic_TB)["assay_raw"]
leic_reproc <- assays(leic_TB)["assay_reprocess"]
leic_reproc <- leic_reproc$assay_reprocess
leic_meta <- colData(leic_TB)
dge_leic <- DGEList(leic_reproc, group = leic_meta$TBStatus)
keep <- rowSums(cpm(dge_leic) > 2) >= 5
summary(keep)
#   Mode   FALSE    TRUE 
#logical  12163   13206 
dge_leic <- dge_leic[keep, , keep.lib.sizes = F]
cpm_leic <- cpm(dge_leic, log = F)
write.csv(cpm_leic, "leicester_cpms.csv")
lcpm_leic <- cpm(dge_leic, log = F)
write.csv(lcpm_leic, "leicester_lcpms.csv")
dge_leic <- calcNormFactors(dge_leic, method = "TMM")
des_leic <- model.matrix(~ dge_leic$samples$group)
dge_leic <- estimateDisp(dge_leic, design = des_leic)
fit_leic <- glmFit(dge_leic, des_leic)
lrt.leic_latent <- glmLRT(fit_leic, coef = 2)
lrt.leic_PTB <- glmLRT(fit_leic, coef = 3)
res.leic_latent <- topTags(lrt.leic_latent, n = nrow(lrt.leic_latent), adjust.method = "BH", sort.by = "PValue")
res.leic_latent <- as.data.frame(res.leic_latent)
leic_latent_genes <- rownames(res.leic_latent)
lgs <- ensembl[match(leic_latent_genes, ensembl$gene_name),]
res.leic_latent$gene_id <- lgs$gene_id
res.leic_latent$transcript_id <- lgs$transcript_id
res.leic_latent$entrez <- lgs$entrez
res.leic_latent$gene_biotype <- lgs$gene_biotype
res.leic_latent$chromosome <- lgs$chromosome
write.csv(res.leic_latent, "full_res_leic_latent.csv")
res.leic_PTB <- topTags(lrt.leic_PTB, n = nrow(lrt.leic_PTB), adjust.method = "BH", sort.by = "PValue")
res.leic_PTB <- as.data.frame(res.leic_PTB)
leic_PTB_genes <- rownames(res.leic_PTB)
lgs <- ensembl[match(leic_PTB_genes, ensembl$gene_name),]
res.leic_PTB$gene_id <- lgs$gene_id
res.leic_PTB$transcript_id <- lgs$transcript_id
res.leic_PTB$entrez <- lgs$entrez
res.leic_PTB$gene_biotype <- lgs$gene_biotype
res.leic_PTB$chromosome <- lgs$chromosome
write.csv(res.leic_PTB, "full_res_leic_PTB.csv")
sig.res.leic_latent <- as.data.frame(res.leic_latent)
sig.res.leic_latent <- sig.res.leic_latent[sig.res.leic_latent$FDR < 0.05,]
sig.res.leic_latent <- sig.res.leic_latent[abs(sig.res.leic_latent$logFC) > 0.25,]
write.csv(sig.res.leic_latent, "sig_res_leic_latent.csv")
sig.res.leic_PTB <- as.data.frame(res.leic_PTB)
sig.res.leic_PTB <- sig.res.leic_PTB[sig.res.leic_PTB$FDR < 0.05,]
sig.res.leic_PTB <- sig.res.leic_PTB[abs(sig.res.leic_PTB$logFC) > 0.25,]
write.csv(sig.res.leic_PTB, "sig_res_leic_PTB.csv")

library(Hmisc)
wilk_comb <- read.csv("wilk_avg_lfc_by_participant.csv", header = T)
colnames(wilk_comb) <- c("gene_name","logFC")
wilk_comb_max <- summarize(wilk_comb$logFC, wilk_comb$gene_name, max)
colnames(wilk_comb_max) <- c("gene_name","logFC_max")
wilk_comb_min <- summarize(wilk_comb$logFC, wilk_comb$gene_name, min)
colnames(wilk_comb_min) <- c("gene_name","logFC_min")
wilk_comb_both <- wilk_comb_max
wilk_comb_both$logFC_min <- wilk_comb_min$logFC_min
wilk_comb_both$logFC <- wilk_comb_max$logFC_max
wilk_comb_both <- as.data.frame(wilk_comb_both)
for(i in 1:nrow(wilk_comb_both)){
  if(abs(wilk_comb_both[i,2]) > abs(wilk_comb_both[i,3])){
    wilk_comb_both[i,4] <- wilk_comb_both[i,2]
    } else if(abs(wilk_comb_both[i,2]) < abs(wilk_comb_both[i,3])){
      wilk_comb_both[i,4] <- wilk_comb_both[i,3]
    }
}
wilk_comb_both <- wilk_comb_both[,c("gene_name","logFC")]
rownames(wilk_comb_both) <- wilk_comb_both$gene_name
write.csv(wilk_comb_both, "wilk_comb.csv")
rownames(wilk_comb_both) <- c()

sig.res.leic_latent$gene_name <- rownames(sig.res.leic_latent)
sig.res.leic_latent <- as.data.frame(sig.res.leic_latent[,c("gene_name","logFC")])
rownames(sig.res.leic_latent) <- c()
wilk_comp <- wilk_comb_both[wilk_comb_both$gene_name %in% sig.res.leic_latent$gene_name,]
leic_latent_comp <- sig.res.leic_latent[sig.res.leic_latent$gene_name %in% wilk_comb_both$gene_name,]
wilk_v_leic_latent_comp <- merge(wilk_comp, leic_latent_comp, by = "gene_name")

cor.test(wilk_v_leic_latent_comp$logFC.x, wilk_v_leic_latent_comp$logFC.y, method = "spearman")

sig.res.leic_PTB$gene_name <- rownames(sig.res.leic_PTB)
sig.res.leic_PTB <- as.data.frame(sig.res.leic_PTB[,c("gene_name","logFC")])
rownames(sig.res.leic_PTB) <- c()
wilk_comp <- wilk_comb_both[wilk_comb_both$gene_name %in% sig.res.leic_PTB$gene_name,]
leic_PTB_comp <- sig.res.leic_PTB[sig.res.leic_PTB$gene_name %in% wilk_comb_both$gene_name,]
wilk_v_leic_PTB_comp <- merge(wilk_comp, leic_PTB_comp, by = "gene_name")

cor.test(wilk_v_leic_PTB_comp$logFC.x, wilk_v_leic_PTB_comp$logFC.y, method = "spearman")

## Make a log counts, CPM and log CPM assay
hivtb_data <- TB_hiv
hivtb_data <- mkAssay(hivtb_data, log = TRUE, counts_to_CPM = TRUE)

### Check to see that we now have 4 assays
assays(hivtb_data)

siglist_hivtb <- names(COVID_sigs)
out_hivtb <- capture.output(ssgsea_result <- runTBsigProfiler(input = hivtb_data,
                                                        useAssay = "log_cpm",
                                                        signatures = COVID_sigs,
                                                        algorithm = "ssGSEA",
                                                        combineSigAndAlgorithm = TRUE,
                                                        parallel.sz = 1))

colors <- RColorBrewer::brewer.pal(6, "Spectral")
col.me <- circlize::colorRamp2(seq(from = -2, to = 2,
                                   length.out = 6), colors)

signatureHeatmap(ssgsea_result, name = "Heatmap of Signatures, 
                 ssGSEA Algorithm",
                 signatureColNames = names(COVID_sigs),
                 annotationColNames = "Disease",
                 scale = TRUE,
                 showColumnNames = TRUE,
                 choose_color = col.me)

signatureBoxplot(inputData = ssgsea_result,
                 name = "Boxplots of Signatures, ssGSEA",
                 signatureColNames = names(COVID_sigs),
                 annotationColName = "Disease", rotateLabels = FALSE)

tb <- read.csv("ProcessedDataMatrix_counts.csv", header = T, row.names = 1)
gene_id <- rownames(tb)
my_genes <- ensembl[match(gene_id, ensembl$gene_id),]
my_genes <- my_genes[!duplicated(my_genes$gene_name),]
my_genes <- my_genes[!is.na(my_genes$gene_name),]
gene_id <- my_genes$gene_id
tb <- tb[gene_id,]
rownames(tb) <- my_genes$gene_name
tb_meta <- read.table("meta.txt", sep = "\t", header = T)
tb_meta <- unique(tb_meta)
tb_meta <- tb_meta[-361,]

tb_SE <- SummarizedExperiment(assays = list(counts = as.matrix(tb)), colData = tb_meta)
assays(tb_SE)

tb_scores_GSVA <- runTBsigProfiler(tb_SE, useAssay = "counts", signatures = COVID_sigs, algorithm = "GSVA")

signatureBoxplot(inputData = tb_scores_GSVA, 
                 name = "GSVA",
                 signatureColNames = names(COVID_sigs),
                 annotationColName = "group", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

tb_SE <- mkAssay(tc_SE, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)

wilk_comb <- read.csv("wilk_comb_max.csv", row.names = 1, header = T)
wilk_genes <- rownames(wilk_comb)
my_genes <- ensembl[match(wilk_genes, ensembl$gene_name),]
my_genes <- my_genes[!duplicated(my_genes$gene_name),]
my_genes <- my_genes[!is.na(my_genes$gene_name),]
gene_names <- my_genes$gene_name
wilk_comb <- as.data.frame(wilk_comb[gene_names,])
rownames(wilk_comb) <- my_genes$gene_id
write.csv(wilk_comb, "wilk_comb_idb.csv")

#gene ontology analyses
#plot InnateDB results
library(ggplot2)
library(VennDiagram)
gg_go_ora <- function(x){
  ggplot(x, aes(x= log.P.value, y= Pathway.Name)) +
    geom_point(aes(colour = rev(log.P.value), size= as.numeric(Pathway.uploaded.gene.count))) +
    scale_color_gradientn(colours=rainbow(4),limits=c(2,15)) +
    geom_vline(xintercept=1.30103, size=0.5, colour="gray50") +
    theme(panel.background = element_rect(fill="gray95", colour="gray95"),
          axis.title.y=element_blank()) +
    scale_y_discrete(limits=rev(reorder(x$Pathway.Name, x$log.P.value)))
}

go_leic_latent <- read.csv("idb_gene_ont_leic_latent.csv", header = T)[,1:8]
go_leic_latent$Source.Name <- as.character(go_leic_latent$Source.Name)
go_leic_latent$Pathway.p.value <- as.numeric(as.character(go_leic_latent$Pathway.p.value))
go_leic_latent_sig <- go_leic_latent[go_leic_latent$Pathway.p.value..corrected.<0.05,]
go_leic_latent_bp <- go_leic_latent[grep("biological process",go_leic_latent$Source.Name),]
go_leic_latent_bp <- go_leic_latent_bp[go_leic_latent_bp$Pathway.p.value..corrected.<0.05,]
go_leic_latent_bp$log.P.value <- -log10(go_leic_latent_bp$Pathway.p.value)
go_leic_latent_mf <- go_leic_latent[grep("molecular function",go_leic_latent$Source.Name),]
go_leic_latent_mf <- go_leic_latent_mf[go_leic_latent_mf$Pathway.p.value..corrected.<0.05,]
go_leic_latent_mf$log.P.value <- -log10(go_leic_latent_mf$Pathway.p.value)
go_leic_latent_cc <- go_leic_latent[grep("cellular component",go_leic_latent$Source.Name),]
go_leic_latent_cc <- go_leic_latent_cc[go_leic_latent_cc$Pathway.p.value..corrected.<0.05,]
go_leic_latent_cc$log.P.value <- -log10(go_leic_latent_cc$Pathway.p.value)

gg_leic_latent_bp <- gg_go_ora(go_leic_latent_bp)
gg_leic_latent_mf <- gg_go_ora(go_leic_latent_mf)
gg_leic_latent_cc <- gg_go_ora(go_leic_latent_cc)

go_leic_PTB <- read.csv("idb_gene_ont_leic_PTB.csv", header = T)[,1:8]
go_leic_PTB$Source.Name <- as.character(go_leic_PTB$Source.Name)
go_leic_PTB$Pathway.p.value <- as.numeric(as.character(go_leic_PTB$Pathway.p.value))
go_leic_PTB_sig <- go_leic_PTB[go_leic_PTB$Pathway.p.value..corrected.<0.05,]
go_leic_PTB_bp <- go_leic_PTB[grep("biological process",go_leic_PTB$Source.Name),]
go_leic_PTB_bp <- go_leic_PTB_bp[go_leic_PTB_bp$Pathway.p.value..corrected.<0.05,]
go_leic_PTB_bp$log.P.value <- -log10(go_leic_PTB_bp$Pathway.p.value)
go_leic_PTB_mf <- go_leic_PTB[grep("molecular function",go_leic_PTB$Source.Name),]
go_leic_PTB_mf <- go_leic_PTB_mf[go_leic_PTB_mf$Pathway.p.value..corrected.<0.05,]
go_leic_PTB_mf$log.P.value <- -log10(go_leic_PTB_mf$Pathway.p.value)
go_leic_PTB_cc <- go_leic_PTB[grep("cellular component",go_leic_PTB$Source.Name),]
go_leic_PTB_cc <- go_leic_PTB_cc[go_leic_PTB_cc$Pathway.p.value..corrected.<0.05,]
go_leic_PTB_cc$log.P.value <- -log10(go_leic_PTB_cc$Pathway.p.value)

gg_leic_PTB_bp <- gg_go_ora(go_leic_PTB_bp)
gg_leic_PTB_mf <- gg_go_ora(go_leic_PTB_mf)
gg_leic_PTB_cc <- gg_go_ora(go_leic_PTB_cc)

go_wilk_comb <- read.csv("idb_gene_ont_wilk_comb.csv", header = T)[,1:8]
go_wilk_comb$Source.Name <- as.character(go_wilk_comb$Source.Name)
go_wilk_comb$Pathway.p.value <- as.numeric(as.character(go_wilk_comb$Pathway.p.value))
go_wilk_comb_sig <- go_wilk_comb[go_wilk_comb$Pathway.p.value..corrected.<0.01,]
go_wilk_comb_bp <- go_wilk_comb[grep("biological process",go_wilk_comb$Source.Name),]
go_wilk_comb_bp <- go_wilk_comb_bp[go_wilk_comb_bp$Pathway.p.value..corrected.<0.001,]
go_wilk_comb_bp$log.P.value <- -log10(go_wilk_comb_bp$Pathway.p.value)
go_wilk_comb_mf <- go_wilk_comb[grep("molecular function",go_wilk_comb$Source.Name),]
go_wilk_comb_mf <- go_wilk_comb_mf[go_wilk_comb_mf$Pathway.p.value..corrected.<0.01,]
go_wilk_comb_mf$log.P.value <- -log10(go_wilk_comb_mf$Pathway.p.value)
go_wilk_comb_cc <- go_wilk_comb[grep("cellular component",go_wilk_comb$Source.Name),]
go_wilk_comb_cc <- go_wilk_comb_cc[go_wilk_comb_cc$Pathway.p.value..corrected.<0.01,]
go_wilk_comb_cc$log.P.value <- -log10(go_wilk_comb_cc$Pathway.p.value)

gg_wilk_comb_bp <- gg_go_ora(go_wilk_comb_bp)
gg_wilk_comb_mf <- gg_go_ora(go_wilk_comb_mf)
gg_wilk_comb_cc <- gg_go_ora(go_wilk_comb_cc)

go_sig_overlap <- venn(list(go_leic_latent_sig$Pathway.Name, go_wilk_comb_sig$Pathway.Name, go_leic_PTB_sig$Pathway.Name))
grid.newpage()                                        
draw.triple.venn(area1 = 14,                          
                 area2 = 87,
                 area3 = 213,
                 n12 = 7,
                 n23 = 55,
                 n13 = 9,
                 n123 = 6,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))
go_bp_overlap <- venn(list(go_leic_latent_bp$Pathway.Name, go_wilk_comb_bp$Pathway.Name, go_leic_PTB_bp$Pathway.Name))
grid.newpage()
draw.triple.venn(area1 = 9,                          
                 area2 = 60,
                 area3 = 126,
                 n12 = 4,
                 n23 = 38,
                 n13 = 6,
                 n123 = 4,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))
go_mf_overlap <- venn(list(go_leic_latent_mf$Pathway.Name, go_wilk_comb_mf$Pathway.Name, go_leic_PTB_mf$Pathway.Name))
grid.newpage()
draw.triple.venn(area1 = 1,                          
                 area2 = 10,
                 area3 = 38,
                 n12 = 1,
                 n23 = 4,
                 n13 = 0,
                 n123 = 0,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))
go_cc_overlap <- venn(list(go_leic_latent_cc$Pathway.Name, go_wilk_comb_cc$Pathway.Name, go_leic_PTB_cc$Pathway.Name))
grid.newpage()
draw.triple.venn(area1 = 4,                          
                 area2 = 17,
                 area3 = 49,
                 n12 = 2,
                 n23 = 13,
                 n13 = 3,
                 n123 = 2,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))

path_leic_latent <- read.csv("idb_pathway_leic_latent.csv", header = T)[,1:8]
path_leic_latent$Source.Name <- as.character(path_leic_latent$Source.Name)
path_leic_latent$Pathway.p.value <- as.numeric(as.character(path_leic_latent$Pathway.p.value))
path_leic_latent_sig <- path_leic_latent[path_leic_latent$Pathway.p.value..corrected.<0.01,]
path_leic_latent_react <- path_leic_latent[grep("REACTOME",path_leic_latent$Source.Name),]
path_leic_latent_react <- path_leic_latent_react[path_leic_latent_react$Pathway.p.value..corrected.<0.01,]
path_leic_latent_react$log.P.value <- -log10(path_leic_latent_react$Pathway.p.value)
path_leic_latent_inoh <- path_leic_latent[grep("INOH",path_leic_latent$Source.Name),]
path_leic_latent_inoh <- path_leic_latent_inoh[path_leic_latent_inoh$Pathway.p.value..corrected.<0.1,]
path_leic_latent_inoh$log.P.value <- -log10(path_leic_latent_inoh$Pathway.p.value)
path_leic_latent_kegg <- path_leic_latent[grep("KEGG",path_leic_latent$Source.Name),]
path_leic_latent_kegg <- path_leic_latent_kegg[path_leic_latent_kegg$Pathway.p.value..corrected.<0.05,]
path_leic_latent_kegg$log.P.value <- -log10(path_leic_latent_kegg$Pathway.p.value)

path_leic_PTB <- read.csv("idb_pathway_leic_PTB.csv", header = T)[,1:8]
path_leic_PTB$Source.Name <- as.character(path_leic_PTB$Source.Name)
path_leic_PTB$Pathway.p.value <- as.numeric(as.character(path_leic_PTB$Pathway.p.value))
path_leic_PTB_sig <- path_leic_PTB[path_leic_PTB$Pathway.p.value..corrected.<0.05,]
path_leic_PTB_react <- path_leic_PTB[grep("REACTOME",path_leic_PTB$Source.Name),]
path_leic_PTB_react <- path_leic_PTB_react[path_leic_PTB_react$Pathway.p.value..corrected.<0.05,]
path_leic_PTB_react$log.P.value <- -log10(path_leic_PTB_react$Pathway.p.value)
path_leic_PTB_inoh <- path_leic_PTB[grep("INOH",path_leic_PTB$Source.Name),]
path_leic_PTB_inoh <- path_leic_PTB_inoh[path_leic_PTB_inoh$Pathway.p.value..corrected.<0.1,]
path_leic_PTB_inoh$log.P.value <- -log10(path_leic_PTB_inoh$Pathway.p.value)
path_leic_PTB_kegg <- path_leic_PTB[grep("KEGG",path_leic_PTB$Source.Name),]
path_leic_PTB_kegg <- path_leic_PTB_kegg[path_leic_PTB_kegg$Pathway.p.value..corrected.<0.05,]
path_leic_PTB_kegg$log.P.value <- -log10(path_leic_PTB_kegg$Pathway.p.value)

path_wilk_comb <- read.csv("idb_pathway_wilk_comb.csv", header = T)[,1:8]
path_wilk_comb$Source.Name <- as.character(path_wilk_comb$Source.Name)
path_wilk_comb$Pathway.p.value <- as.numeric(as.character(path_wilk_comb$Pathway.p.value))
path_wilk_comb_sig <- path_wilk_comb[path_wilk_comb$Pathway.p.value..corrected.<0.001,]
path_wilk_comb_react <- path_wilk_comb[grep("REACTOME",path_wilk_comb$Source.Name),]
path_wilk_comb_react <- path_wilk_comb_react[path_wilk_comb_react$Pathway.p.value..corrected.<0.001,]
path_wilk_comb_react$log.P.value <- -log10(path_wilk_comb_react$Pathway.p.value)
path_wilk_comb_inoh <- path_wilk_comb[grep("INOH",path_wilk_comb$Source.Name),]
path_wilk_comb_inoh <- path_wilk_comb_inoh[path_wilk_comb_inoh$Pathway.p.value..corrected.<0.1,]
path_wilk_comb_inoh$log.P.value <- -log10(path_wilk_comb_inoh$Pathway.p.value)
path_wilk_comb_kegg <- path_wilk_comb[grep("KEGG",path_wilk_comb$Source.Name),]
path_wilk_comb_kegg <- path_wilk_comb_kegg[path_wilk_comb_kegg$Pathway.p.value..corrected.<0.01,]
path_wilk_comb_kegg$log.P.value <- -log10(path_wilk_comb_kegg$Pathway.p.value)

gg_react_ora <- function(x){
  ggplot(data = x, aes(x=reorder(Pathway.Name, log.P.value), y= log.P.value)) +
    geom_bar(stat = "identity") +
    geom_col(aes(fill = log.P.value)) +
    scale_fill_gradient2(low = "firebrick1", high = "firebrick4", name = "-log10(P-value)") +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.5) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
}
gg_inoh_ora <- function(x){
  ggplot(data = x, aes(x=reorder(Pathway.Name, log.P.value), y= log.P.value)) +
    geom_bar(stat = "identity") +
    geom_col(aes(fill = log.P.value)) +
    scale_fill_gradient2(low = "dodgerblue", high = "dodgerblue4", name = "-log10(P-value)") +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.5) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
}
gg_kegg_ora <- function(x){
  ggplot(data = x, aes(x=reorder(Pathway.Name, log.P.value), y= log.P.value)) +
    geom_bar(stat = "identity") +
    geom_col(aes(fill = log.P.value)) +
    scale_fill_gradient2(low = "darkolivegreen1", high = "darkolivegreen4", name = "-log10(P-value)") +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.5) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
}

gg_path_leic_latent_react <- gg_react_ora(path_leic_latent_react)
gg_path_leic_latent_inoh <- gg_inoh_ora(path_leic_latent_inoh)
gg_path_leic_latent_kegg <- gg_kegg_ora(path_leic_latent_kegg)

gg_path_leic_PTB_react <- gg_react_ora(path_leic_PTB_react)
gg_path_leic_PTB_inoh <- gg_inoh_ora(path_leic_PTB_inoh)
gg_path_leic_PTB_kegg <- gg_kegg_ora(path_leic_PTB_kegg)

gg_path_wilk_comb_react <- gg_react_ora(path_wilk_comb_react)
gg_path_wilk_comb_inoh <- gg_inoh_ora(path_wilk_comb_inoh)
gg_path_wilk_comb_kegg <- gg_kegg_ora(path_wilk_comb_kegg)

path_sig_overlap <- venn(list(path_leic_latent_sig$Pathway.Name, path_wilk_comb_sig$Pathway.Name, path_leic_PTB_sig$Pathway.Name))
grid.newpage()
draw.triple.venn(area1 = 24,                          
                 area2 = 54,
                 area3 = 89,
                 n12 = 9,
                 n23 = 29,
                 n13 = 13,
                 n123 = 8,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))
path_react_overlap <- venn(list(path_leic_latent_react$Pathway.Name, path_wilk_comb_react$Pathway.Name, path_leic_PTB_react$Pathway.Name))
grid.newpage()
draw.triple.venn(area1 = 17,                          
                 area2 = 25,
                 area3 = 42,
                 n12 = 7,
                 n23 = 14,
                 n13 = 11,
                 n123 = 6,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))
path_inoh_overlap <- venn(list(path_leic_latent_inoh$Pathway.Name, path_wilk_comb_inoh$Pathway.Name, path_leic_PTB_inoh$Pathway.Name))
grid.newpage()
draw.triple.venn(area1 = 2,                          
                 area2 = 2,
                 area3 = 6,
                 n12 = 0,
                 n23 = 1,
                 n13 = 0,
                 n123 = 0,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))
path_kegg_overlap <- venn(list(path_leic_latent_kegg$Pathway.Name, path_wilk_comb_kegg$Pathway.Name, path_leic_PTB_kegg$Pathway.Name))
grid.newpage()
draw.triple.venn(area1 = 9,                          
                 area2 = 12,
                 area3 = 20,
                 n12 = 4,
                 n23 = 8,
                 n13 = 5,
                 n123 = 3,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))

tf_leic_latent <- read.csv("idb_TFBS_leic_latent.csv", header = T)[,1:8]
tf_leic_latent$Source.Name <- as.character(tf_leic_latent$Source.Name)
tf_leic_latent$Pathway.p.value <- as.numeric(as.character(tf_leic_latent$Pathway.p.value))
tf_leic_latent_sig <- tf_leic_latent[tf_leic_latent$Pathway.p.value<0.05,]
tf_leic_latent_sig$log.P.value <- -log10(tf_leic_latent_sig$Pathway.p.value)

tf_leic_PTB <- read.csv("idb_TFBS_leic_PTB.csv", header = T)[,1:8]
tf_leic_PTB$Source.Name <- as.character(tf_leic_PTB$Source.Name)
tf_leic_PTB$Pathway.p.value <- as.numeric(as.character(tf_leic_PTB$Pathway.p.value))
tf_leic_PTB_sig <- tf_leic_PTB[tf_leic_PTB$Pathway.p.value<0.05,]
tf_leic_PTB_sig$log.P.value <- -log10(tf_leic_PTB_sig$Pathway.p.value)

tf_wilk_comb <- read.csv("idb_TFBS_wilk_comb.csv", header = T)[,1:8]
tf_wilk_comb$Source.Name <- as.character(tf_wilk_comb$Source.Name)
tf_wilk_comb$Pathway.p.value <- as.numeric(as.character(tf_wilk_comb$Pathway.p.value))
tf_wilk_comb_sig <- tf_wilk_comb[tf_wilk_comb$Pathway.p.value<0.05,]
tf_wilk_comb_sig$log.P.value <- -log10(tf_wilk_comb_sig$Pathway.p.value)

gg_tfbs_ora <- function(x){
  ggplot(data = x, aes(x=reorder(Pathway.Name, log.P.value), y= log.P.value)) +
    geom_bar(stat = "identity") +
    geom_col(aes(fill = log.P.value)) +
    scale_fill_gradient2(low = "mediumpurple1", high = "mediumpurple4", name = "-log10(P-value)") +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.5) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
}

gg_tf_leic_latent <- gg_tfbs_ora(tf_leic_latent_sig)
gg_tf_leic_PTB <- gg_tfbs_ora(tf_leic_PTB_sig)
gg_tf_wilk_comb <- gg_tfbs_ora(tf_wilk_comb_sig)

tf_overlap <- venn(list(tf_leic_latent_sig$Pathway.Name, tf_wilk_comb_sig$Pathway.Name, tf_leic_PTB_sig$Pathway.Name))
grid.newpage()
draw.triple.venn(area1 = 9,                          
                 area2 = 20,
                 area3 = 29,
                 n12 = 4,
                 n23 = 7,
                 n13 = 2,
                 n123 = 2,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Latent TB", "Pulmonary TB", "COVID-19"))


#prepare O'Garra influenza dataset
library(GEOquery)
gset <- getGEO("GSE111368", GSEMatrix =TRUE)
targets <- pData(gset[[1]])
setwd("/Users/sheerin.d/OneDrive\ -\ wehi.edu.au/COVID-19/meta_analysis/GSE111368")
x <- read.ilmn("/Users/sheerin.d/OneDrive\ -\ wehi.edu.au/COVID-19/meta_analysis/GSE111368/GSE111368_Non-normalized_data.txt",probeid="ID_REF")
y <- neqc(x)
colnames(y) <- rownames(targets)

#Build the design matrix for the linear modelling function
f <- factor(targets$`flu_type:ch1`, levels = unique(targets$`flu_type:ch1`))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
rownames(design) <- rownames(targets)

#Apply the intensity values to lmFit
fit <- lmFit(y, design)
write.table(fit, file="fit.txt", sep="\t", quote=FALSE)

#Create a contrast matrix
contrast.matrix <- makeContrasts("A-HC", levels=design)

#Apply this contrast matrix to the modeled data and compute statistics for the data
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Output the statistics for the dataset and write them to disk
output <- topTreat(fit2, coef=1, number=Inf, adjust.method="BH", lfc=0.58)
output <- output[output$adj.P.Val < 0.05,]
ez_id <- rownames(output)
my_ez <- ensembl[match(ez_id, ensembl$entrez),]
my_ez <- my_ez[!duplicated(my_ez$gene_name),]
output <- output[rownames(output)%in%my_ez$entrez,]
my_genes <- my_ez$gene_name[-1]
rownames(output) <- my_genes
ilmn <- rownames(y)
ilmn <- output[match(ilmn, output$ID),]
ilmn <- ilmn[!duplicated(rownames(ilmn)),]
ilmn <- ilmn[!is.na(rownames(ilmn)),]
y <- y[rownames(y) %in% output$ID,]
ilmn <- ilmn[ilmn$ID %in% rownames(y),]
rownames(y) <- rownames(ilmn)
write.table(output[,c("logFC","adj.P.Val")], file="ogarra_flu.csv", sep=",", quote=FALSE)

flu_SE <- SummarizedExperiment(assays = list(lcpm = as.matrix(y)), colData = targets)

#prepare Silvin COVID-19 whole blood data
silvin <- read.csv("silvin_DEGs.csv")
silvin$p_val_adj <- as.numeric(silvin$p_val_adj)
silvin <- silvin[silvin$p_val_adj < 0.05,]
silvin <- silvin[abs(silvin$avg_logFC) > 0.58,]
write.csv(silvin, "silvin.csv")

#get O'Garra TB data
ogarra_SE <- readRDS(file = "/Users/sheerin.d/OneDrive - wehi.edu.au/2020/Bioinformatics/HHC_study/kristina_analysis/ogarra_SE.rds")
ogarra_tb <- assays(ogarra_SE)["counts"]$counts
ogarra_meta <- as.data.frame(colData(ogarra_SE))
tb_status <- ogarra_meta$TBStatus
tb_status[tb_status == "Active_TB"] = "Leicester_Active"
tb_status[tb_status == "Leicester_LTBI"] = "Leicester_Latent"
tb_status <- factor(tb_status)
tb_status <- relevel(tb_status, ref = "Leicester_Control")
dge_og <- DGEList(ogarra_tb, group = tb_status)
keep <- rowSums(cpm(dge_og) > 2) >= 5
summary(keep)
#   Mode   FALSE    TRUE 
#logical   43734   16908
dge_og <- dge_og[keep, , keep.lib.sizes = F]
cpm_og <- cpm(dge_og, log = F)
write.csv(cpm_og, "ogarra_TB_cpms.csv")
lcpm_og <- cpm(dge_og, log = F)
write.csv(lcpm_og, "ogarra_TB_lcpms.csv")
dge_og <- calcNormFactors(dge_og, method = "TMM")
des_og <- model.matrix(~ + dge_og$samples$group)
dge_og <- estimateDisp(dge_og, design = des_og)
fit_og <- glmFit(dge_og, des_og)
lrt.og_leic_active <- glmLRT(fit_og, coef = 2)
lrt.og_leic_latent <- glmLRT(fit_og, coef = 3)
lrt.og_leic_prog <- glmLRT(fit_og, coef = 4)
lrt.og_lon_latent <- glmLRT(fit_og, coef = 5)
lrt.og_lon_PTB <- glmLRT(fit_og, coef = 6)
lrt.og_SA_latent <- glmLRT(fit_og, coef = 7)
lrt.og_SA_PTB <- glmLRT(fit_og, coef = 7)
res.og_leic_active <- topTags(lrt.og_leic_active, n = nrow(lrt.og_leic_active), adjust.method = "BH", sort.by = "PValue")
res.og_leic_active <- as.data.frame(res.og_leic_active)
og_leic_active_genes <- rownames(res.og_leic_active)
lgs <- ensembl[match(og_leic_active_genes, ensembl$gene_name),]
res.og_leic_active$gene_id <- lgs$gene_id
res.og_leic_active$transcript_id <- lgs$transcript_id
res.og_leic_active$entrez <- lgs$entrez
res.og_leic_active$gene_biotype <- lgs$gene_biotype
res.og_leic_active$chromosome <- lgs$chromosome
write.csv(res.og_leic_active, "full_res_og_leic_active.csv")
res.og_leic_latent <- topTags(lrt.og_leic_latent, n = nrow(lrt.og_leic_latent), adjust.method = "BH", sort.by = "PValue")
res.og_leic_latent <- as.data.frame(res.og_leic_latent)
og_leic_latent_genes <- rownames(res.og_leic_latent)
lgs <- ensembl[match(og_leic_latent_genes, ensembl$gene_name),]
res.og_leic_latent$gene_id <- lgs$gene_id
res.og_leic_latent$transcript_id <- lgs$transcript_id
res.og_leic_latent$entrez <- lgs$entrez
res.og_leic_latent$gene_biotype <- lgs$gene_biotype
res.og_leic_latent$chromosome <- lgs$chromosome
write.csv(res.og_leic_latent, "full_res_og_leic_latent.csv")
res.og_leic_prog <- topTags(lrt.og_leic_prog, n = nrow(lrt.og_leic_prog), adjust.method = "BH", sort.by = "PValue")
res.og_leic_prog <- as.data.frame(res.og_leic_prog)
og_leic_prog_genes <- rownames(res.og_leic_prog)
lgs <- ensembl[match(og_leic_prog_genes, ensembl$gene_name),]
res.og_leic_prog$gene_id <- lgs$gene_id
res.og_leic_prog$transcript_id <- lgs$transcript_id
res.og_leic_prog$entrez <- lgs$entrez
res.og_leic_prog$gene_biotype <- lgs$gene_biotype
res.og_leic_prog$chromosome <- lgs$chromosome
write.csv(res.og_leic_prog, "full_res_og_leic_prog.csv")
res.og_lon_latent <- topTags(lrt.og_lon_latent, n = nrow(lrt.og_lon_latent), adjust.method = "BH", sort.by = "PValue")
res.og_lon_latent <- as.data.frame(res.og_lon_latent)
og_lon_latent_genes <- rownames(res.og_lon_latent)
lgs <- ensembl[match(og_lon_latent_genes, ensembl$gene_name),]
res.og_lon_latent$gene_id <- lgs$gene_id
res.og_lon_latent$transcript_id <- lgs$transcript_id
res.og_lon_latent$entrez <- lgs$entrez
res.og_lon_latent$gene_biotype <- lgs$gene_biotype
res.og_lon_latent$chromosome <- lgs$chromosome
write.csv(res.og_lon_latent, "full_res_og_lon_latent.csv")
res.og_lon_PTB <- topTags(lrt.og_lon_PTB, n = nrow(lrt.og_lon_PTB), adjust.method = "BH", sort.by = "PValue")
res.og_lon_PTB <- as.data.frame(res.og_lon_PTB)
og_lon_PTB_genes <- rownames(res.og_lon_PTB)
lgs <- ensembl[match(og_lon_PTB_genes, ensembl$gene_name),]
res.og_lon_PTB$gene_id <- lgs$gene_id
res.og_lon_PTB$transcript_id <- lgs$transcript_id
res.og_lon_PTB$entrez <- lgs$entrez
res.og_lon_PTB$gene_biotype <- lgs$gene_biotype
res.og_lon_PTB$chromosome <- lgs$chromosome
write.csv(res.og_lon_PTB, "full_res_og_lon_PTB.csv")
res.og_SA_latent <- topTags(lrt.og_SA_latent, n = nrow(lrt.og_SA_latent), adjust.method = "BH", sort.by = "PValue")
res.og_SA_latent <- as.data.frame(res.og_SA_latent)
og_SA_latent_genes <- rownames(res.og_SA_latent)
lgs <- ensembl[match(og_SA_latent_genes, ensembl$gene_name),]
res.og_SA_latent$gene_id <- lgs$gene_id
res.og_SA_latent$transcript_id <- lgs$transcript_id
res.og_SA_latent$entrez <- lgs$entrez
res.og_SA_latent$gene_biotype <- lgs$gene_biotype
res.og_SA_latent$chromosome <- lgs$chromosome
write.csv(res.og_SA_latent, "full_res_og_SA_latent.csv")
res.og_SA_PTB <- topTags(lrt.og_SA_PTB, n = nrow(lrt.og_SA_PTB), adjust.method = "BH", sort.by = "PValue")
res.og_SA_PTB <- as.data.frame(res.og_SA_PTB)
og_SA_PTB_genes <- rownames(res.og_SA_PTB)
lgs <- ensembl[match(og_SA_PTB_genes, ensembl$gene_name),]
res.og_SA_PTB$gene_id <- lgs$gene_id
res.og_SA_PTB$transcript_id <- lgs$transcript_id
res.og_SA_PTB$entrez <- lgs$entrez
res.og_SA_PTB$gene_biotype <- lgs$gene_biotype
res.og_SA_PTB$chromosome <- lgs$chromosome
write.csv(res.og_SA_PTB, "full_res_og_SA_PTB.csv")

sig.res.og_leic_active <- as.data.frame(res.og_leic_active)
sig.res.og_leic_active <- sig.res.og_leic_active[sig.res.og_leic_active$FDR < 0.05,]
sig.res.og_leic_active <- sig.res.og_leic_active[abs(sig.res.og_leic_active$logFC) > 0.58,]
write.csv(sig.res.og_leic_active[,c("logFC", "FDR")], "sig_res_og_leic_active.csv")
sig.res.og_leic_latent <- as.data.frame(res.og_leic_latent)
sig.res.og_leic_latent <- sig.res.og_leic_latent[sig.res.og_leic_latent$FDR < 0.05,]
sig.res.og_leic_latent <- sig.res.og_leic_latent[abs(sig.res.og_leic_latent$logFC) > 0.25,]
write.csv(sig.res.og_leic_latent[,c("logFC", "FDR")], "sig_res_og_leic_latent.csv")
sig.res.og_leic_prog <- as.data.frame(res.og_leic_prog)
sig.res.og_leic_prog <- sig.res.og_leic_prog[sig.res.og_leic_prog$FDR < 0.05,]
sig.res.og_leic_prog <- sig.res.og_leic_prog[abs(sig.res.og_leic_prog$logFC) > 0.25,]
write.csv(sig.res.og_leic_prog[,c("logFC", "FDR")], "sig_res_og_leic_prog.csv")
sig.res.og_lon_latent <- as.data.frame(res.og_lon_latent)
sig.res.og_lon_latent <- sig.res.og_leic_latent[sig.res.og_lon_latent$FDR < 0.05,]
sig.res.og_lon_latent <- sig.res.og_lon_latent[abs(sig.res.og_lon_latent$logFC) > 0.25,]
write.csv(sig.res.og_lon_latent[,c("logFC", "FDR")], "sig_res_og_lon_latent.csv")
sig.res.og_lon_PTB <- as.data.frame(res.og_lon_PTB)
sig.res.og_lon_PTB <- sig.res.og_leic_prog[sig.res.og_lon_PTB$FDR < 0.05,]
sig.res.og_lon_PTB <- sig.res.og_leic_prog[abs(sig.res.og_lon_PTB$logFC) > 0.25,]
write.csv(sig.res.og_lon_PTB[,c("logFC", "FDR")], "sig_res_og_lon_PTB.csv")
sig.res.og_SA_latent <- as.data.frame(res.og_SA_latent)
sig.res.og_SA_latent <- sig.res.og_SA_latent[sig.res.og_SA_latent$FDR < 0.05,]
sig.res.og_SA_latent <- sig.res.og_SA_latent[abs(sig.res.og_SA_latent$logFC) > 0.25,]
write.csv(sig.res.og_SA_latent[,c("logFC", "FDR")], "sig_res_og_SA_latent.csv")
sig.res.og_SA_PTB <- as.data.frame(res.og_SA_PTB)
sig.res.og_SA_PTB <- sig.res.og_SA_PTB[sig.res.og_SA_PTB$FDR < 0.05,]
sig.res.og_SA_PTB <- sig.res.og_SA_PTB[abs(sig.res.og_SA_PTB$logFC) > 0.25,]
write.csv(sig.res.og_SA_PTB[,c("logFC", "FDR")], "sig_res_og_SA_PTB.csv")


#get Pulendran COVID-19 PBMC bulk RNA-seq data
bali_rna <- getGEO("GSE152418", GSEMatrix = TRUE)
bali_meta <- pData(bali_rna[[1]])
bali_meta <- bali_meta[,c(1,44:50)]
colnames(bali_meta) <- c("title","supp_file","age","cell_type","disease_state","sex","location","severity")
rownames(bali_meta) <- bali_meta$title
bali_meta <- bali_meta[,-1]
bali <- read.table("/Users/sheerin.d/OneDrive - wehi.edu.au/2020/COVID-19/meta_analysis/Pulendran_PBMCs_bulk/GSE152418_p20047_Study1_RawCounts.txt", header = T, row.names = 1)
rownames(bali_meta) <- gsub("-",".",rownames(bali_meta))
bali_meta <- bali_meta[-1,]
bali <- bali[,-1]

all(rownames(bali_meta) %in% colnames(bali))
#TRUE
all(rownames(bali_meta) == colnames(bali))
#TRUE
disease_state <- factor(bali_meta$disease_state)
disease_state <- relevel(disease_state, ref = "Healthy")
severity <- factor(bali_meta$severity)
severity <- relevel(severity, ref = "Healthy")

#create summarizedExperiment for Bali data
bali_ids <- rownames(bali)
bali_ids <- ensembl[match(bali_ids, ensembl$gene_id),]
bali_ids <- bali_ids[!duplicated(bali_ids$gene_name),]
bali <- bali[rownames(bali) %in% bali_ids$gene_id,]
bali_ids <- bali_ids[bali_ids$gene_id %in% rownames(bali),]
rownames(bali) <- bali_ids$gene_name

bali_SE <- SummarizedExperiment(assays = list(counts = as.matrix(bali)), colData = bali_meta)
bali_SE <- mkAssay(bali_SE, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)

dge_bali <- DGEList(bali, group = disease_state)
keep <- rowSums(cpm(dge_bali) > 0) >= 17
summary(keep)
#   Mode   FALSE    TRUE 
#logical   34209   26474
dge_bali <- dge_bali[keep, , keep.lib.sizes = F]
dge_bali <- calcNormFactors(dge_bali, method = "TMM")
des_bali <- model.matrix(~ + dge_bali$samples$group)
dge_bali <- estimateDisp(dge_bali, design = des_bali)
fit_bali <- glmFit(dge_bali, des_bali)
lrt.cov_v_hea <- glmLRT(fit_bali, coef = 2)
res.cov_v_hea <- topTags(lrt.cov_v_hea, n = nrow(lrt.cov_v_hea), adjust.method = "BH", sort.by = "PValue")
res.cov_v_hea <- as.data.frame(res.cov_v_hea)
cov_v_hea_genes <- rownames(res.cov_v_hea)
gs <- ensembl[match(cov_v_hea_genes, ensembl$gene_id),]
gs <- gs[!duplicated(gs$gene_name),]
res.cov_v_hea <- res.cov_v_hea[rownames(res.cov_v_hea)%in%gs$gene_id,]
rownames(res.cov_v_hea) <- gs$gene_name
sig.res.cov_v_hea <- res.cov_v_hea[res.cov_v_hea$FDR < 0.01,]
sig.res.cov_v_hea <- sig.res.cov_v_hea[abs(sig.res.cov_v_hea$logFC) > 0.58,]
write.csv(sig.res.cov_v_hea[,c("logFC", "FDR")], "sig_res_COVID_v_healthy_BALI_bulk_PBMCs.csv")

des_bali <- model.matrix(~ + severity)
dge_bali <- estimateDisp(dge_bali, design = des_bali)
fit_bali <- glmFit(dge_bali, des_bali)
lrt.icu <- glmLRT(fit_bali, coef = 2)
lrt.mod <- glmLRT(fit_bali, coef = 3)
lrt.sev <- glmLRT(fit_bali, coef = 4)

res.icu <- topTags(lrt.icu, n = nrow(lrt.icu), adjust.method = "BH", sort.by = "PValue")
res.mod <- topTags(lrt.mod, n = nrow(lrt.mod), adjust.method = "BH", sort.by = "PValue")
res.sev <- topTags(lrt.sev, n = nrow(lrt.sev), adjust.method = "BH", sort.by = "PValue")
res.icu <- as.data.frame(res.icu)
res.mod <- as.data.frame(res.mod)
res.sev <- as.data.frame(res.sev)

icu_genes <- rownames(res.icu)
gs <- ensembl[match(icu_genes, ensembl$gene_id),]
gs <- gs[!duplicated(gs$gene_name),]
res.icu <- res.icu[rownames(res.icu)%in%gs$gene_id,]
rownames(res.icu) <- gs$gene_name

mod_genes <- rownames(res.mod)
gs <- ensembl[match(mod_genes, ensembl$gene_id),]
gs <- gs[!duplicated(gs$gene_name),]
res.mod <- res.mod[rownames(res.mod)%in%gs$gene_id,]
rownames(res.mod) <- gs$gene_name

sev_genes <- rownames(res.sev)
gs <- ensembl[match(sev_genes, ensembl$gene_id),]
gs <- gs[!duplicated(gs$gene_name),]
res.sev <- res.sev[rownames(res.sev)%in%gs$gene_id,]
rownames(res.sev) <- gs$gene_name

sig.res.icu <- res.icu[res.icu$FDR < 0.01,]
sig.res.icu <- sig.res.icu[abs(sig.res.icu$logFC) > 0.58,]
write.csv(sig.res.icu[,c("logFC", "FDR")], "sig_res_ICU_BALI_bulk_PBMCs.csv")

sig.res.mod <- res.mod[res.mod$FDR < 0.01,]
sig.res.mod <- sig.res.mod[abs(sig.res.mod$logFC) > 0.58,]
write.csv(sig.res.mod[,c("logFC", "FDR")], "sig_res_MOD_BALI_bulk_PBMCs.csv")

sig.res.sev <- res.sev[res.sev$FDR < 0.01,]
sig.res.sev <- sig.res.sev[abs(sig.res.sev$logFC) > 0.58,]
write.csv(sig.res.sev[,c("logFC", "FDR")], "sig_res_SEV_BALI_bulk_PBMCs.csv")

#Bali using COVID-19/TB overlapping pathway genes
COVID_spec1 <- read.csv("Pathway_99.csv")
COVID_spec2 <- read.csv("Silvin_511_genes.csv")
COVID_spec3 <- read.csv("20_99_selection.csv")
COVID_spec4 <- read.csv("65_Arunachalam_enrichment_Silvin.csv")
COVID_spec5 <- read.csv("only_COVID.csv")
COVID_spec <- list(as.character(COVID_spec1[,1]), as.character(COVID_spec2[,1]), as.character(COVID_spec3[,1]), as.character(COVID_spec4[,1]), as.character(COVID_spec5[,1]))
names(COVID_spec) <- c("Pathway_99","Silvin_511","C20_99_selection","C65_Arunachalam_enrichment_Silvin","only_COVID")

bali_scores_GSVA <- runTBsigProfiler(bali_SE, useAssay = "log_cpm", signatures = COVID_spec, algorithm = "GSVA")
og_scores_GSVA <- runTBsigProfiler(ogarra_SE, useAssay = "log_cpm", signatures = COVID_spec, algorithm = "GSVA")

signatureBoxplot(inputData = bali_scores_GSVA, 
                 name = "GSVA",
                 signatureColNames = names(COVID_spec),
                 annotationColName = "severity", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

signatureBoxplot(inputData = og_scores_GSVA, 
                 name = "GSVA",
                 signatureColNames = names(COVID_spec),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


#Silvin DE pathways
deg_silvin <- read.csv("DEG_list_SIlvin_pathways.csv", header = T)
covid_atb_ltb_no_flu <- as.character(deg_silvin$COVID_ATB_LTB_NO_FLU)[1:93]
covid_atb_flu_no_ltb <- as.character(deg_silvin$COVID_ATB_FLU_NO_LTB)[1:213]
covid_atb_no_ltb_no_flu <- as.character(deg_silvin$COVID_ATB_NO_LTB_NO_FLU)
deg_silvin <- list(COVID_ATB_LTB_no_flu = covid_atb_flu_no_ltb, COVID_ATB_FLU_no_ltb = covid_atb_ltb_no_flu, COVID_ATB_no_ltb_no_flu = covid_atb_no_ltb_no_flu)
silvin_mtorc <- read.csv("MTORC_signaling_Silvin.csv")
silvin_mtorc <- as.character(silvin_mtorc[,1])
silvin <- list(silvin_mtorc = silvin_mtorc, silvin_20 = as.character(COVID_spec3[,1]))

gsva_PTB_Latent <- lapply(multi_set_PTB_Latent,
                           function(x) TBSignatureProfiler::runTBsigProfiler(
                             input = x,
                             useAssay = assayNames(x),
                             signatures = deg_silvin,
                             algorithm = "GSVA",
                             combineSigAndAlgorithm = TRUE))

gsva_PTB_Latent1 <- lapply(multi_set_PTB_Latent,
                          function(x) TBSignatureProfiler::runTBsigProfiler(
                            input = x,
                            useAssay = assayNames(x),
                            signatures = silvin,
                            algorithm = "GSVA",
                            combineSigAndAlgorithm = TRUE))

bali_scores_GSVA1 <- runTBsigProfiler(bali_SE, useAssay = "log_cpm", signatures = silvin, algorithm = "GSVA")

flu_scores_GSVA <- runTBsigProfiler(flu_SE, useAssay = "lcpm", signatures = silvin_mtorc, algorithm = "GSVA")

ogarra_scores_GSVA1 <- runTBsigProfiler(ogarra_SE, useAssay = "log_cpm", signatures = silvin, algorithm = "GSVA")

signatureBoxplot(inputData = ogarra_scores_GSVA1, 
                 name = "GSVA",
                 signatureColNames = names(silvin),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#look at all the Leicester groups first
colnames(colData(gsva_PTB_Latent1$GSE107994))
both_GSVA <- as.data.frame(colData(gsva_PTB_Latent1$GSE107994)[,c(4,27,28)])
colnames(both_GSVA) <- c("Condition","silvin_mtorc","silvin_20")
colnames(colData(bali_scores_GSVA1))
bali_scores_GSVA <- as.data.frame(colData(bali_scores_GSVA1)[,c(7:9)])
bali_scores_GSVA$severity <- factor(bali_scores_GSVA$severity)
levels(bali_scores_GSVA$severity[bali_scores_GSVA$severity == "Healthy"]) = "1. Healthy (Arunachalam et al.)"
bali_scores_GSVA$severity <- factor(bali_scores_GSVA$severity, levels = c("1. Healthy (Arunachalam et al.)","10. Moderate COVID-19","11. Severe COVID-19","ICU COVID-19"))
colnames(bali_scores_GSVA) <- c("Condition","silvin_mtorc","silvin_20")
both_GSVA <- rbind(both_GSVA, bali_scores_GSVA)
both_GSVA$Condition <- factor(both_GSVA$Condition, levels = c("Latent","PTB","Healthy","Moderate","Severe","ICU"))
both_GSVA <- reshape2::melt(both_GSVA, id.vars = c("Condition"))
levels(both_GSVA$variable) <- c("COVID-19 mTORC Signalling Pathway","Select COVID-19 Pathway Genes")
ggplot(both_GSVA, aes(x = Condition, y = value, group = Condition)) +
  geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~ variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Latent", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14), legend.position = "bottom")

both_GSVA <- as.data.frame(colData(ogarra_scores_GSVA1)[,c(4,23,24)])
colnames(both_GSVA) <- c("Condition","silvin_mtorc","silvin_20")
colnames(colData(bali_scores_GSVA1))
bali_scores_GSVA <- as.data.frame(colData(bali_scores_GSVA1)[,c(7:9)])
bali_scores_GSVA$severity <- factor(bali_scores_GSVA$severity, levels = c("9. Healthy","10. Moderate","11. Severe","12. ICU"))
colnames(bali_scores_GSVA) <- c("Condition","silvin_mtorc","silvin_20")
both_GSVA <- rbind(both_GSVA, bali_scores_GSVA)
both_GSVA$Condition <- factor(both_GSVA$Condition, levels = c("Latent","PTB","Healthy","Moderate","Severe","ICU"))
both_GSVA <- reshape2::melt(both_GSVA, id.vars = c("Condition"))
levels(both_GSVA$variable) <- c("COVID-19 mTORC Signalling Pathway","Select COVID-19 Pathway Genes")
ggplot(both_GSVA, aes(x = Condition, y = value, group = Condition)) +
  geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~ variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Healthy", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14), legend.position = "bottom")

bali_scores_GSVA <- runTBsigProfiler(bali_SE, useAssay = "log_cpm", signatures = deg_silvin, algorithm = "GSVA")
og_scores_GSVA <- runTBsigProfiler(ogarra_SE, useAssay = "log_cpm", signatures = deg_silvin, algorithm = "GSVA")

colnames(colData(og_scores_GSVA))
og_scores_GSVA <- as.data.frame(colData(og_scores_GSVA)[,c(4,23:25)])
og_scores_GSVA <- reshape2::melt(og_scores_GSVA, id.vars = c("TBStatus"))
levels(og_scores_GSVA$variable) <- c("Present in COVID-19, Active and Latent TB (absent in Influenza)","Present in COVID-19, Active TB and Influenza (absent in Latent TB)","Present in COVID-19 and Active TB (absent in Latent TB and Influenza")
ggplot(og_scores_GSVA, aes(x = TBStatus, y = value, group = TBStatus)) +
  geom_boxplot(aes(fill = TBStatus), outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~ variable, nrow = 3) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "1. London Latent", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14), legend.position = "bottom")

stat.test01 <- both_GSVA %>%
  group_by(variable) %>%
  t_test(value ~ Condition, ref.group = "Latent") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test01 <- as.data.frame(stat.test01)
write.csv(stat.test01, "GSEA_lists_OG_stats.csv")

stat.test01 <- both_GSVA %>%
  group_by(variable) %>%
  t_test(value ~ Condition, ref.group = "Latent") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test01 <- as.data.frame(stat.test01)
write.csv(stat.test01, "GSEA_lists_OG_stats.csv")

signatureBoxplot(inputData = bali_scores_GSVA, 
                 name = "Gene Set Variation Analysis of Arunachalam et. al. COVID-19 data",
                 signatureColNames = names(deg_silvin),
                 annotationColName = "severity", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
signatureBoxplot(inputData = og_scores_GSVA, 
                 name = "Gene Set Variation Analysis of O'Garra et. al. TB data",
                 signatureColNames = names(deg_silvin),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
signatureBoxplot(inputData = og_scores_GSVA, 
                 name = "Gene Set Variation Analysis of O'Garra et. al. TB data",
                 signatureColNames = names(silvin_mtorc),
                 annotationColName = "TBStatus", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
bali_scores_GSVA <- runTBsigProfiler(bali_SE, useAssay = "log_cpm", signatures = silvin_mtorc, algorithm = "GSVA")
signatureBoxplot(inputData = bali_scores_GSVA, 
                 name = "Gene Set Variation Analysis of Pulendran et. al. COVID-19 data",
                 signatureColNames = names(silvin_mtorc),
                 annotationColName = "severity", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

bali_scores_GSVA <- runTBsigProfiler(bali_SE, useAssay = "log_cpm", signatures = silvin_mtorc, algorithm = "GSVA")
og_scores_GSVA <- runTBsigProfiler(ogarra_SE, useAssay = "log_cpm", signatures = silvin_mtorc, algorithm = "GSVA")

colnames(colData(og_scores_GSVA))
og_scores_GSVA <- as.data.frame(colData(og_scores_GSVA)[,c(4,23)])
og_scores_GSVA <- reshape2::melt(og_scores_GSVA, id.vars = c("TBStatus"))
levels(og_scores_GSVA$variable) <- c("COVID-19 mTORC1 Signalling Pathway (TB)")
ggplot(og_scores_GSVA, aes(x = TBStatus, y = value, group = TBStatus)) +
  geom_boxplot(aes(fill = TBStatus), outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~ variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "1. London Latent", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

stat.test02 <- og_scores_GSVA %>%
  group_by(variable) %>%
  t_test(value ~ TBStatus, ref.group = "1. London Latent") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test02 <- as.data.frame(stat.test02)
write.csv(stat.test02, "mTORC1_OG_stats.csv")

colnames(colData(bali_scores_GSVA))
bali_scores_GSVA <- as.data.frame(colData(bali_scores_GSVA)[,c(7,8)])
bali_scores_GSVA <- reshape2::melt(bali_scores_GSVA, id.vars = c("severity"))
levels(bali_scores_GSVA$variable) <- c("COVID-19 mTORC1 Signalling Pathway (Arunachalam et. al.)")
bali_scores_GSVA$severity <- factor(bali_scores_GSVA$severity, levels = c("Healthy","Moderate","Severe","ICU"))
ggplot(bali_scores_GSVA, aes(x = severity, y = value, group = severity)) +
  geom_boxplot(aes(fill = severity), outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~ variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Healthy", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

stat.test03 <- bali_scores_GSVA %>%
  group_by(variable) %>%
  t_test(value ~ severity, ref.group = "Healthy") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test03 <- as.data.frame(stat.test03)
write.csv(stat.test03, "mTORC1_bali_stats.csv")

#HHC data
fc_base_SE <- SummarizedExperiment(assays = list(counts = as.matrix(fc_base)), colData = meta_base)
fc_base_SE <- mkAssay(fc_base_SE, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)

#run scoring algorithm
HHC_scores_GSVA <- runTBsigProfiler(fc_base_SE, useAssay = "log_cpm", signatures = COVID_rearr, algorithm = "GSVA")
HHC_GSVA_res_grp <- as.data.frame(colData(HHC_scores_GSVA)[,c(7,29:64)])
HHC_GSVA_res_par <- as.data.frame(colData(HHC_scores_GSVA)[,c(16,29:64)])
HHC_GSVA_res_par$par[HHC_GSVA_res_par$par == "no_lesion"] = "1. No lesion"
HHC_GSVA_res_par$par[HHC_GSVA_res_par$par == "cold"] = "2. Other cold par"
HHC_GSVA_res_par$par[HHC_GSVA_res_par$par == "def1"] = "3. Subclinical definition 1"
HHC_GSVA_res_par$par[HHC_GSVA_res_par$par == "def2"] = "4. Subclinical definition 2"
HHC_GSVA_res_par$par[HHC_GSVA_res_par$par == "hot"] = "5. Other hot par"
HHC_GSVA_res_ln_lung <- as.data.frame(colData(HHC_scores_GSVA)[,c(20,29:64)])
HHC_GSVA_res_ln_lung$ln_lung[HHC_GSVA_res_ln_lung$ln_lung == "neither"] = "1. Neither"
HHC_GSVA_res_ln_lung$ln_lung[HHC_GSVA_res_ln_lung$ln_lung == "ln"] = "2. Hot lymph node"
HHC_GSVA_res_ln_lung$ln_lung[HHC_GSVA_res_ln_lung$ln_lung == "lung"] = "3. Hot lung"
HHC_GSVA_res_ln_lung$ln_lung[HHC_GSVA_res_ln_lung$ln_lung == "both"] = "4. Both"
HHC_GSVA_res_ch <- as.data.frame(colData(HHC_scores_GSVA)[,c(21,29:64)])
HHC_GSVA_res_ch <- HHC_GSVA_res_ch[grep("no change|improvement|worsening",HHC_GSVA_res_ch$change),]
HHC_GSVA_res_ch$change[HHC_GSVA_res_ch$change == "no change"] = "1. No change"
HHC_GSVA_res_ch$change[HHC_GSVA_res_ch$change == "improvement"] = "2. Improvement"
HHC_GSVA_res_ch$change[HHC_GSVA_res_ch$change == "worsening"] = "3. Worsening"
HHC_GSVA_res_stat <- as.data.frame(colData(HHC_scores_GSVA)[,c(22,29:64)])
HHC_GSVA_res_stat <- HHC_GSVA_res_stat[grep("no TB|incident|prevalent",HHC_GSVA_res_stat$status),]
HHC_GSVA_res_stat$status[HHC_GSVA_res_stat$status == "no TB"] = "1. No TB"
HHC_GSVA_res_stat$status[HHC_GSVA_res_stat$status == "incident"] = "2. Incident"
HHC_GSVA_res_stat$status[HHC_GSVA_res_stat$status == "prevalent"] = "3. Prevalent"

HHC_GSVA_res_grp <- reshape2::melt(HHC_GSVA_res_grp, id.vars = c("group"))
levels(HHC_GSVA_res_grp$variable) <- c("Wilk Monocytes", "Wen Monocytes", "Wen Inactivated Monocytes", "Wei Classical Monocytes", "Silvin Monocytes (WB)", "Wilk NK Cells", "Liao G1 M  (BALF)", "Liao G1/2 M  (BALF)", "Liao G2 M  (BALF)", "Liao G3 M  (BALF)", "Liao G4 M  (BALF)", "Wen NK cells", "Wen CD4 T cells", "Wen CD8 T cells", "Liao CD8 T cells", "Wei T cells", "Wen B Cells", "Wei B Cells", "Wilk Activated Granulocytes", "Silvin Neutrophils (WB)", "Silvin WB Signature", "Wilk ISG Signature", "Huang IFN Signature", "Hadjadj ISG Signature (WB)", "Hadjadj WB Signature", "Hadjadj Mild/Moderate (WB)", "Hadjadj Severe (WB)", "Hadjadj Critical (WB)", "Xiong PBMC Signature", "Xiong Cytokine Signature", "Arunachalam PBMC Signature", "Arunachalam Moderate", "Arunachalam Severe", "Arunachalam ICU", "Arunachalam BTMs", "Dunning Influenza (WB)")
HHC_GSVA_res_par <- reshape2::melt(HHC_GSVA_res_par, id.vars = c("par"))
levels(HHC_GSVA_res_par$variable) <- c("Wilk Monocytes", "Wen Monocytes", "Wen Inactivated Monocytes", "Wei Classical Monocytes", "Silvin Monocytes (WB)", "Wilk NK Cells", "Liao G1 M  (BALF)", "Liao G1/2 M  (BALF)", "Liao G2 M  (BALF)", "Liao G3 M  (BALF)", "Liao G4 M  (BALF)", "Wen NK cells", "Wen CD4 T cells", "Wen CD8 T cells", "Liao CD8 T cells", "Wei T cells", "Wen B Cells", "Wei B Cells", "Wilk Activated Granulocytes", "Silvin Neutrophils (WB)", "Silvin WB Signature", "Wilk ISG Signature", "Huang IFN Signature", "Hadjadj ISG Signature (WB)", "Hadjadj WB Signature", "Hadjadj Mild/Moderate (WB)", "Hadjadj Severe (WB)", "Hadjadj Critical (WB)", "Xiong PBMC Signature", "Xiong Cytokine Signature", "Arunachalam PBMC Signature", "Arunachalam Moderate", "Arunachalam Severe", "Arunachalam ICU", "Arunachalam BTMs", "Dunning Influenza (WB)")
HHC_GSVA_res_ln_lung <- reshape2::melt(HHC_GSVA_res_ln_lung, id.vars = c("ln_lung"))
levels(HHC_GSVA_res_ln_lung$variable) <- c("Wilk Monocytes", "Wen Monocytes", "Wen Inactivated Monocytes", "Wei Classical Monocytes", "Silvin Monocytes (WB)", "Wilk NK Cells", "Liao G1 M  (BALF)", "Liao G1/2 M  (BALF)", "Liao G2 M  (BALF)", "Liao G3 M  (BALF)", "Liao G4 M  (BALF)", "Wen NK cells", "Wen CD4 T cells", "Wen CD8 T cells", "Liao CD8 T cells", "Wei T cells", "Wen B Cells", "Wei B Cells", "Wilk Activated Granulocytes", "Silvin Neutrophils (WB)", "Silvin WB Signature", "Wilk ISG Signature", "Huang IFN Signature", "Hadjadj ISG Signature (WB)", "Hadjadj WB Signature", "Hadjadj Mild/Moderate (WB)", "Hadjadj Severe (WB)", "Hadjadj Critical (WB)", "Xiong PBMC Signature", "Xiong Cytokine Signature", "Arunachalam PBMC Signature", "Arunachalam Moderate", "Arunachalam Severe", "Arunachalam ICU", "Arunachalam BTMs", "Dunning Influenza (WB)")
HHC_GSVA_res_ch <- reshape2::melt(HHC_GSVA_res_ch, id.vars = c("change"))
levels(HHC_GSVA_res_ch$variable) <- c("Wilk Monocytes", "Wen Monocytes", "Wen Inactivated Monocytes", "Wei Classical Monocytes", "Silvin Monocytes (WB)", "Wilk NK Cells", "Liao G1 M  (BALF)", "Liao G1/2 M  (BALF)", "Liao G2 M  (BALF)", "Liao G3 M  (BALF)", "Liao G4 M  (BALF)", "Wen NK cells", "Wen CD4 T cells", "Wen CD8 T cells", "Liao CD8 T cells", "Wei T cells", "Wen B Cells", "Wei B Cells", "Wilk Activated Granulocytes", "Silvin Neutrophils (WB)", "Silvin WB Signature", "Wilk ISG Signature", "Huang IFN Signature", "Hadjadj ISG Signature (WB)", "Hadjadj WB Signature", "Hadjadj Mild/Moderate (WB)", "Hadjadj Severe (WB)", "Hadjadj Critical (WB)", "Xiong PBMC Signature", "Xiong Cytokine Signature", "Arunachalam PBMC Signature", "Arunachalam Moderate", "Arunachalam Severe", "Arunachalam ICU", "Arunachalam BTMs", "Dunning Influenza (WB)")
HHC_GSVA_res_stat <- reshape2::melt(HHC_GSVA_res_stat, id.vars = c("status"))
levels(HHC_GSVA_res_stat$variable) <- c("Wilk Monocytes", "Wen Monocytes", "Wen Inactivated Monocytes", "Wei Classical Monocytes", "Silvin Monocytes (WB)", "Wilk NK Cells", "Liao G1 M  (BALF)", "Liao G1/2 M  (BALF)", "Liao G2 M  (BALF)", "Liao G3 M  (BALF)", "Liao G4 M  (BALF)", "Wen NK cells", "Wen CD4 T cells", "Wen CD8 T cells", "Liao CD8 T cells", "Wei T cells", "Wen B Cells", "Wei B Cells", "Wilk Activated Granulocytes", "Silvin Neutrophils (WB)", "Silvin WB Signature", "Wilk ISG Signature", "Huang IFN Signature", "Hadjadj ISG Signature (WB)", "Hadjadj WB Signature", "Hadjadj Mild/Moderate (WB)", "Hadjadj Severe (WB)", "Hadjadj Critical (WB)", "Xiong PBMC Signature", "Xiong Cytokine Signature", "Arunachalam PBMC Signature", "Arunachalam Moderate", "Arunachalam Severe", "Arunachalam ICU", "Arunachalam BTMs", "Dunning Influenza (WB)")

#group
ggplot(HHC_GSVA_res_grp, aes(x = group, y = value, group = group)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA) +
  geom_jitter(aes(), width = 0.2, size = 0.5) +
  facet_wrap(~variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  scale_fill_discrete(name = "Group") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "1", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

#par
ggplot(HHC_GSVA_res_par, aes(x = par, y = value, group = par)) +
  geom_boxplot(aes(fill = par), outlier.shape = NA) +
  geom_jitter(aes(), width = 0.2, size = 0.5) +
  facet_wrap(~variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  scale_fill_discrete(name = "Group") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "1. No lesion", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

#LN/lung
ggplot(HHC_GSVA_res_ln_lung, aes(x = ln_lung, y = value, group = ln_lung)) +
  geom_boxplot(aes(fill = ln_lung), outlier.shape = NA) +
  geom_jitter(aes(), width = 0.2, size = 0.5) +
  facet_wrap(~variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  scale_fill_discrete(name = "Group") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "1. Neither", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

#change
ggplot(HHC_GSVA_res_ch, aes(x = change, y = value, group = change)) +
  geom_boxplot(aes(fill = change), outlier.shape = NA) +
  geom_jitter(aes(), width = 0.2, size = 0.5) +
  facet_wrap(~variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  scale_fill_discrete(name = "Group") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "1. No change", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

#status
ggplot(HHC_GSVA_res_stat, aes(x = status, y = value, group = status)) +
  geom_boxplot(aes(fill = status), outlier.shape = NA) +
  geom_jitter(aes(), width = 0.2, size = 0.5) +
  facet_wrap(~variable) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  scale_fill_discrete(name = "Group") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "1. No TB", hide.ns = T, label.y = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

