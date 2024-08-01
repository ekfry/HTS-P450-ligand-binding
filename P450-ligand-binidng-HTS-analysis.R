
#     The following script was written to process spectral data obtained from a plate reader as referred to in 
#  "Development of a high throughput cytochrome P450-ligand binding assay" published in XX. The script is equipped to 
#  analyze cytochrome P450 binding assay data from high throughput screens. The script corrects spectra for background
#  absorbance, prints absolute and difference spectra for each assay, and estimated dissociation constants and ∆Amax
#  values for each assay based on the one-site, specific binding model.




#install.packages("tidyverse")
#install.packages("devtools")
#devtools::install_github("padpadpadpad/nls.multstart")
#install.packages("nls.multstart")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(nls.multstart)
library(broom)


############# (step 1) export spectral information from the plate reader as a csv file #############

############# (step 2) read raw data file, normalize all spectra to start from same absorbance value at 500 nm, and print corrected spectra #############

# insert name of raw data file
r = read.csv("2D6_day94_titrations_RT.csv")

#define protein concentration
protein_conc <- 1
rr = r %>% select(PlateNumber, Well, Wavelength, Result)
r2 =  rr %>% mutate(absorbance = rr$Result) %>% select(-"Result")
r2$wellid <- paste(r2$Plate,r2$Well)
r2$wellid = as.factor(r2$wellid)
# r = r %>% select("wavelength", "wellid", "absorbance")
r3 = reshape(r2, idvar = "Wavelength", timevar = "wellid", v.names = "absorbance", direction = "wide")
names(r3) <- make.names(names(r3), unique=TRUE)
r4 = select(r3, c(-1,-2,-3))
anchor = list()
for (x in 2:ncol(r4)) {
  anchor[x-1] = r4[301, x] - r4[301,1]
}
for (x in 2:ncol(r4)) {
  for (i in 1:301) {
    r4[i, x] = r4[i, x] - anchor[[x-1]]
  }
}

write_csv(r4, "anchored.csv")

############# (step 2) measure background control wells (C1, C2, P1, and P2) and water-bound P450 spectra (A1, A2, N1, and N2) for each plate in the assay  #############

baseline1 <- rowMeans(r4[, c("absorbance.1.C01", "absorbance.1.C02", "absorbance.1.P01", "absorbance.1.P02")])
abs1 <- rowMeans(r4[, c("absorbance.1.A01", "absorbance.1.A02", "absorbance.1.N01", "absorbance.1.N02")])
abs1_cor <- abs1 - baseline1
baseline2 <- rowMeans(r4[, c("absorbance.2.C01", "absorbance.2.C02", "absorbance.2.P01", "absorbance.2.P02")])
abs2 <- rowMeans(r4[, c("absorbance.2.A01", "absorbance.2.A02", "absorbance.2.N01", "absorbance.2.N02")])
abs2_cor <- abs2 - baseline2

#subtract average background absorbance as measured by wells C1, C2, P1, and P2 from spectra in the corresponding plate
baseline_cor1 = r4[, c(1:107)] - baseline1
baseline_cor2 = r4[, c(108:467)] - baseline2

#subtract average water-bound P450 absorbance as measured by wells A1, A2, N1, and N2 from spectra in the corresponding plate
baseline_cor1_d = baseline_cor1[, c(1:107)] - abs1_cor
baseline_cor2_d = baseline_cor2[, c(1:360)] - abs2_cor

#print all corrected absolute spectra
abs <- cbind(baseline_cor2, baseline_cor3, baseline_cor4) %>% mutate(wavelength = seq(350,500,0.5)) %>% relocate(wavelength)
write_csv(abs, "abs.csv")

#print all corrected difference spectra
dif <- cbind(baseline_cor2_d, baseline_cor3_d, baseline_cor4_d) %>% mutate(wavelength = seq(350,500,0.5)) %>% relocate(wavelength)
write_csv(dif, "dif.csv")

#edit script to account for as many plates as are present in the assay


############# (step 2) print corrected water-bound P450 spectra for each plate to visually ensure quality of data ############# 

#qc plate 1
baseline1_qc <- data.frame(baseline1)
baseline1_qc1 <- baseline1_qc %>% mutate(wavelength = seq(350,500,0.5)) %>% relocate(wavelength)
abs1_qc_temp <- data.frame(abs1)
abs1_qc <- abs1_qc_temp %>% mutate(wavelength = seq(350,500,0.5)) %>% relocate(wavelength) %>% mutate (plate = 1)
colnames(abs1_qc)[colnames(abs1_qc) == "abs1"] <- "Absorbance"
abs1_cor_qc_temp <- data.frame(abs1_cor)
abs1_cor_qc <- abs1_cor_qc_temp %>% mutate(wavelength = seq(350,500,0.5)) %>% relocate(wavelength) %>% mutate (plate = "1_cor")
colnames(abs1_cor_qc)[colnames(abs1_cor_qc) == "abs1_cor"] <- "Absorbance"
abs1_qc <- rbind(abs1_qc, abs1_cor_qc)

ggplot(baseline1_qc1) + 
  geom_line(aes(x = wavelength, y = baseline1)) + xlim(350, 500) +
  ggtitle("baseline1") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept = 393, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 418, linetype = "dotted") +
  geom_vline(xintercept = 424, linetype = "dotted", color = "red")

ggsave(plot = last_plot(), filename = paste0("baseline1.jpeg"), device = "jpeg", path = "qc",
       width = 6, height = 3, units = "in", dpi = 320)

ggplot(abs1_qc) + 
  geom_line(aes(x = wavelength, y = Absorbance, color = plate, group = plate)) + xlim(350, 500) +
  ggtitle("DMSO control") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept = 393, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 418, linetype = "dotted") +
  geom_vline(xintercept = 424, linetype = "dotted", color = "red")

ggsave(plot = last_plot(), filename = paste0("DMSO_control1.jpeg"), device = "jpeg", path = "qc",
       width = 6, height = 3, units = "in", dpi = 320)

#qc plate 2
dir.create("qc")
baseline2_qc <- data.frame(baseline2)
baseline2_qc2 <- baseline2_qc %>% mutate(wavelength = seq(350,500,0.5)) %>% relocate(wavelength)
abs2_qc_temp <- data.frame(abs2)
abs2_qc <- abs2_qc_temp %>% mutate(wavelength = seq(350,500,0.5)) %>% relocate(wavelength) %>% mutate (plate = 2)
colnames(abs2_qc)[colnames(abs2_qc) == "abs2"] <- "Absorbance"
abs2_cor_qc_temp <- data.frame(abs2_cor)
abs2_cor_qc <- abs2_cor_qc_temp %>% mutate(wavelength = seq(350,500,0.5)) %>% relocate(wavelength) %>% mutate (plate = "2_cor")
colnames(abs2_cor_qc)[colnames(abs2_cor_qc) == "abs2_cor"] <- "Absorbance"
abs2_qc <- rbind(abs2_qc, abs2_cor_qc)

ggplot(baseline2_qc2) + 
  geom_line(aes(x = wavelength, y = baseline2)) + xlim(350, 500) +
  ggtitle("baseline2") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept = 393, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 418, linetype = "dotted") +
  geom_vline(xintercept = 424, linetype = "dotted", color = "red")

ggsave(plot = last_plot(), filename = paste0("baseline2.jpeg"), device = "jpeg", path = "qc",
       width = 6, height = 3, units = "in", dpi = 320)

ggplot(abs2_qc) + 
  geom_line(aes(x = wavelength, y = Absorbance, color = plate, group = plate)) + xlim(350, 500) +
  ggtitle("DMSO control") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept = 393, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 418, linetype = "dotted") +
  geom_vline(xintercept = 424, linetype = "dotted", color = "red")

ggsave(plot = last_plot(), filename = paste0("DMSO_control2.jpeg"), device = "jpeg", path = "qc",
       width = 6, height = 3, units = "in", dpi = 320)

#edit script to account for as many plates as are present in the assay


############# (step 3) generate, group, and print absolute spectra for each assay #############

absolute_spectra <- abs %>% select(-contains("wavelength")) %>% select(-contains("01")) %>% select(-contains("02"))

# define the number of titration points per assay
num_titration_points <- 11

for (i in 1:(ncol(absolute_spectra)/num_titration_points)) {
  temp_index <- i * num_titration_points - (num_titration_points - 1)
  temp_file <- gather(absolute_spectra, key = "ID", value = "Absorbance", c(temp_index:(temp_index + (num_titration_points - 1))))
  temp_file <- temp_file[,c("ID","Absorbance")]
  temp_file$Short.Titration <- i
  
  if (i == 1) {
    final_abs <- temp_file
  } else {
    final_abs <- rbind(final_abs, temp_file)
  }
}

x <- seq(350,500,0.5)
u <- unique(final_abs$Short.Titration)
dir.create("absolute")
abs_DMSO <- rowMeans(abs[, c("absorbance.1.A01", "absorbance.1.A02", "absorbance.1.N01", "absorbance.1.N02")])
abs_DMSO_plot <- data.frame(ID="absorbance.1.DMSO", Absorbance = abs_DMSO, Short.Titration = "x")
compound_key <- read.csv("compound_key_titration_RT.csv")

for (i in 1:length(u)) {
  plot_data <- final_abs %>% filter(Short.Titration == i)
  plot_data1 <- rbind(plot_data, abs_DMSO_plot)
  plot_data2 <- plot_data1 %>% mutate(wavelength = rep(x,(num_titration_points+1)))
  
  ggplot(plot_data2) + 
    geom_line(aes(x = wavelength, y = Absorbance, color = ID, group = ID)) + xlim(350, 500) + ylim(0, 0.1) +
    ggtitle(compound_key[i,1]) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_vline(xintercept = 393, linetype = "dotted", color = "blue") +
    geom_vline(xintercept = 418, linetype = "dotted") +
    geom_vline(xintercept = 424, linetype = "dotted", color = "red")
  
  ggsave(plot = last_plot(), filename = paste0(compound_key[i,1], ".jpeg"), device = "jpeg", path = "absolute",
         width = 6, height = 3, units = "in", dpi = 320)
}


############# (step 3) generate, group, and print difference spectra for each assay #############

difference_spectra <- dif %>% select(-wavelength) %>% select(-contains("01")) %>% select(-contains("02"))

for (i in 1:(ncol(difference_spectra)/num_titration_points)) {
  temp_index <- i * num_titration_points - (num_titration_points - 1)
  temp_file_d <- gather(difference_spectra, key = "ID", value = "Absorbance", c(temp_index:(temp_index+(num_titration_points-1))))
  temp_file_d <- temp_file_d[,c("ID","Absorbance")]
  temp_file_d$Short.Titration <- i
  
  if (i == 1) {
    final_dif <- temp_file_d
  } else {
    final_dif <- rbind(final_dif, temp_file_d)
  }
}

u <- unique(final_dif$Short.Titration)
dir.create("difference")
dif_DMSO <- rowMeans(dif[, c("absorbance.1.A01", "absorbance.1.A02", "absorbance.1.N01", "absorbance.1.N02")])
dif_DMSO_plot <- data.frame(ID="absorbance.1.DMSO", Absorbance = dif_DMSO, Short.Titration = "x")

for (i in 1:length(u)) {
  plot_data_d <- final_dif %>% filter(Short.Titration == i)
  plot_data_d1 <- rbind(plot_data_d, dif_DMSO_plot)
  plot_data_d2 <- plot_data_d1 %>% mutate(wavelength = rep(x,(num_titration_points+1)))
  
  ggplot(plot_data_d2) + 
    geom_line(aes(x = wavelength, y = Absorbance, color = ID, group = ID)) + xlim(350, 500) + ylim(-0.03, 0.03) +
    ggtitle(compound_key[i,1]) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  ggsave(plot = last_plot(), filename = paste0(compound_key[i,1], ".jpeg"), device = "jpeg", path = "difference",
         width = 6, height = 3, units = "in", dpi = 320)
}


############# read in appropriate compound concentrations for each assay ############# 

#read in file containing compound codes or names and concentrations in each well
table_with_conc <- read.csv("240325_2D6 day 60_with_concentrations_RT.csv")
compound_conc <- data.frame(matrix(nrow = 12, ncol = 0))
for (i in 1:nrow(compound_key)) {
  compound_code <- strsplit(compound_key[i,1],"-")[[1]][1]
  # compound codes are the well ID for the master library - ensure no compounds codes loaded for secondary plates
  # must only use each compound code once

  compound_conc[[compound_key[i,1]]] <- c(subset(table_with_conc, `source.well` == compound_code )$`actual..ligand.`, 0)
}


############# (steps 4-5) calculate ∆A and fit each data set to one site, specific binding nonlinear regression #############

dir.create("binding_curves")
dir.create("binding_curves_rawdata")

Kds <- data.frame(Compounds = compound_key, 
                  Kd = NA, 
                  Kd_lci = NA, 
                  Kd_uci = NA, 
                  Amax = NA, 
                  Amax_lci = NA, 
                  Amax_uci = NA,
                  AIC = NA,
                  BIC = NA,
                  peak = NA, 
                  trough = NA)
results <- list()

#group difference spectra for each assay and identify maximum and minimum absorbance values for highest concentration of compound
for (i in 1:(ncol(difference_spectra)/num_titration_points)) {
  temp_index <- i * num_titration_points - (num_titration_points - 1)
  temp_file_t <- difference_spectra[,c(temp_index:(temp_index + (num_titration_points - 1)))]
  temp_file_tt <- temp_file_t %>% mutate(baseline = dif_DMSO)
  temp_file_t1 <- temp_file_tt[62:182,]
  peak <- temp_file_t1[which.max(temp_file_t1[,5]),]
  trough <- temp_file_t1[which.min(temp_file_t1[,5]),]
  dA <- peak - trough
  dA1 <- stack(dA[1:ncol(dA)])
  dA2 <- dA1 %>% select(-ind)
  dA_final <- as.numeric(unlist(dA2))
  
  peak_trough_temp <- temp_file_t1 %>% mutate(wavelength = seq(380,440,0.5))
  peak_nm <- peak_trough_temp$wavelength[which.max(peak_trough_temp[,5])]
  trough_nm <- peak_trough_temp$wavelength[which.min(peak_trough_temp[,5])]
  Kds[i,"peak"] <- peak_nm
  Kds[i,"trough"] <- trough_nm
  
  binding_curve1 <- data.frame(concentration = compound_conc[[compound_key[i,1]]], values = dA_final)
  
  #print the raw values for change in absorbance and concentration of compound for each assay to import into GraphPad Prism to manage data
  write_csv(binding_curve1, file = paste0("binding_curves_rawdata/", compound_key[i,1], ".csv"))
  
  #apply nonlinear regression model one site, specific binding to each assay
  tryCatch({
    one_site_specific_binding <- function(x, Bmax, Kd) {
      y = Bmax * x / (Kd + x)
    }
    
    
    model <- nls_multstart(values ~ one_site_specific_binding(concentration, Bmax, Kd), 
                           data = binding_curve1,
                           iter = 100,
                           start_lower = list(Bmax = max(dA_final), Kd = 0.01),
                           start_upper = list(Bmax = max(dA_final)*1.5, Kd = 25))
    
    
    Amax <- coef(model)["Bmax"]
    Kd <- coef(model)["Kd"]
    Kds[i,"Kd"] <- Kd
    Kds[i,"Amax"] <- Amax
    
    model_statistics <- summary(model)
    model_statistics2 <- glance(model)
    CI <- confint(model)
    Kd_lci <- CI[2,1]
    Kd_uci <- CI[2,2] 
    Amax_lci <- CI[1,1] 
    Amax_uci <- CI[1,2]
    AIC <- model_statistics2$AIC
    BIC <- model_statistics2$BIC
    
    Kds[i,"Kd_lci"] <- Kd_lci
    Kds[i,"Kd_uci"] <- Kd_uci
    Kds[i,"Amax_lci"] <- Amax_lci
    Kds[i,"Amax_uci"] <- Amax_uci
    Kds[i,"AIC"] <- AIC
    Kds[i,"BIC"] <- BIC
    
    curve <- seq(0, max(compound_conc[[compound_key[i,1]]]), length=500)
    
    concentration_test <- data.frame(concentration = curve)
 
    #print each binding curve   
    ggplot() +
      geom_point(aes(x = compound_conc[[compound_key[i,1]]], y = dA_final)) +
      geom_line(aes(x = curve, y = predict(object = model, newdata = concentration_test)), color = "red") +
      xlim(0, max(compound_conc[[compound_key[i,1]]])) +
      annotate("text", x = 5, y = max(dA_final)*1.1, label = paste0("Amax = ", round(Amax, 4), ", Kd = ", round(Kd, 4))) +
      ylim(0, max(dA_final)*1.1) +
      xlab("Concentration") +
      ylab("deltaAbs") +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      ggtitle(compound_key[i,1])
    
    ggsave(plot = last_plot(), filename = paste0(compound_key[i,1], ".jpeg"), device = "jpeg", path = "binding_curves",
           width = 6, height = 3, units = "in", dpi = 320)
    
    write_csv(binding_curve1, file = paste0("binding_curves_rawdata/", compound_key[i,1], ".csv"))
    
    results[[i]] <- model
  }, error = function(err) {
    cat("model fitting failed for data", i, "\n")
  })
  
}

#print a file containing the Kd with upper and lower confidence intervals, Amax with upper and lower confidence intervals, AIC, BIC, and wavelength 
#of the peak and trough absorbance values for each compound
write_csv(Kds, "Kd_Amax_values.csv")


############# (steps 4-5) calculate ∆A and fit each data set to the quadratic formula #############

dir.create("binding_curves_tb")
dir.create("binding_curves_tb_rawdata")
Kds_tb <- data.frame(Compounds = compound_key, 
                  Kd = NA, 
                  Kd_lci = NA, 
                  Kd_uci = NA, 
                  Amax = NA, 
                  Amax_lci = NA, 
                  Amax_uci = NA,
                  AIC = NA,
                  BIC = NA,
                  peak = NA, 
                  trough = NA)
results_tb <- list()


### note: Because the assay was designed with a protein concentration of 1 µM, if the data is 
### consistent with Kd values significantly lower than the protein concentration, one should consider use of the tight binding or 
### Morrison equation

#apply tight binding (Morrison equation) to each assay
for (i in 1:(ncol(difference_spectra)/num_titration_points)) {
  temp_index <- i * num_titration_points - (num_titration_points - 1)
  temp_file_t <- difference_spectra[,c(temp_index:(temp_index + (num_titration_points - 1)))]
  temp_file_tt <- temp_file_t %>% mutate(baseline = dif_DMSO)
  temp_file_t1 <- temp_file_tt[62:182,]
  peak <- temp_file_t1[which.max(temp_file_t1[,5]),]
  trough <- temp_file_t1[which.min(temp_file_t1[,5]),]
  dA <- peak - trough
  dA1 <- stack(dA[1:ncol(dA)])
  dA2 <- dA1 %>% select(-ind)
  dA_final <- as.numeric(unlist(dA2))
  
  peak_trough_temp <- temp_file_t1 %>% mutate(wavelength = seq(380,440,0.5))
  peak_nm <- peak_trough_temp$wavelength[which.max(peak_trough_temp[,5])]
  trough_nm <- peak_trough_temp$wavelength[which.min(peak_trough_temp[,5])]
  Kds_tb[i,"peak"] <- peak_nm
  Kds_tb[i,"trough"] <- trough_nm
  
  binding_curve1 <- data.frame(concentration = compound_conc[[compound_key[i,1]]], values = dA_final)
  
  #apply quadratic formula to each data set
  tryCatch({
    one_site_specific_binding <- function(x, Bmax, Kd) {
      y = (Bmax/2*protein_conc)*((protein_conc + x + Kd) - sqrt((protein_conc + x + Kd)^2 - (4*protein_conc*x)))
    }
    
    
    model <- nls_multstart(values ~ one_site_specific_binding(concentration, Bmax, Kd), 
                           data = binding_curve1,
                           iter = 500,
                           start_lower = c(Bmax = median(dA_final), Kd = 0.01),
                           start_upper = c(Bmax = max(dA_final)*2, Kd = 25),
                           lower = c(Bmax = min(dA_final), Kd = 0))
    
    
    Amax <- coef(model)["Bmax"]
    Kd <- coef(model)["Kd"]
    Kds_tb[i,"Kd"] <- Kd
    Kds_tb[i,"Amax"] <- Amax
    
    model_statistics <- summary(model)
    model_statistics2 <- glance(model)
    CI <- confint(model)
    Kd_lci <- CI[2,1]
    Kd_uci <- CI[2,2] 
    Amax_lci <- CI[1,1] 
    Amax_uci <- CI[1,2]
    AIC <- model_statistics2$AIC
    BIC <- model_statistics2$BIC
    
    Kds_tb[i,"Kd_lci"] <- Kd_lci
    Kds_tb[i,"Kd_uci"] <- Kd_uci
    Kds_tb[i,"Amax_lci"] <- Amax_lci
    Kds_tb[i,"Amax_uci"] <- Amax_uci
    Kds_tb[i,"AIC"] <- AIC
    Kds_tb[i,"BIC"] <- BIC
    
    curve <- seq(0,max(compound_conc[[compound_key[i,1]]]),length=101)
    
    concentration_test <- data.frame(concentration = curve)
    
    #print each binding curve   
    ggplot() +
      geom_point(aes(x = compound_conc[[compound_key[i,1]]], y = dA_final)) +
      geom_line(aes(x = curve, y = predict(object = model, newdata = concentration_test)), color = "red") +
      xlim(0, max(compound_conc[[compound_key[i,1]]])) +
      annotate("text", x = 5, y = max(dA_final)*1.1, label = paste0("Amax = ", round(Amax, 4), ", Kd = ", round(Kd, 4))) +
      ylim(0, max(dA_final)*1.1) +
      xlab("Concentration") +
      ylab("deltaAbs") +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      ggtitle(compound_key[i,1])
    
    ggsave(plot = last_plot(), filename = paste0(compound_key[i,1], ".jpeg"), device = "jpeg", path = "binding_curves_tb",
           width = 6, height = 3, units = "in", dpi = 320)
    write_csv(binding_curve1, file = paste0("binding_curves_tb_rawdata/", compound_key[i,1], ".csv"))
    
    results_tb[[i]] <- model
  }, error = function(err) {
    cat("model fitting failed for data", i, "\n")
  })
  
}

#print a file containing the Kd with upper and lower confidence intervals, Amax with upper and lower confidence intervals, AIC, BIC, and wavelength 
#of the peak and trough absorbance values for each compound
write_csv(Kds_tb, "Kd_Amax_values_tight_binding.csv")


#import data into Prism to check quality
  #under HTS>Prism scripts, scripts have been written to automatically import the data in "binding_curves_rawdata" 
  #and export Kd and Amax values
