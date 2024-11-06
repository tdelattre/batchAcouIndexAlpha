#####################################
#Global acoustic indices calculation#
#####################################

# Authors: amandine gasc: amandine.gasc@gmail.com
# Source material : https://github.com/agasc/Soundscape-analysis-with-R/blob/master/AcouIndexAlpha.r
# adaptations by thomas.delattre@inrae.fr - should be indicated as such below. 
# I'm not taking any credit for anything interesting in here, I just made Amandine's code compatible with batch processing 
# just don't bother Amandine with anything that's been added by me and maybe bugged ;)
# this script is used by calcul_Indices_Acoustiques_batch_SEGMENTS_X.X.R for parallel processing
# please download both code files

#changelog tdelattre : 
## XX.XX.23
# changed how acouIndexAlpha loads inputs to make it compatible with batch processing
## 28.05.24
# modification des fonctions de soundecology pour les rendre silencieuses et accélérer le calcul en batch
# NB : ces fonctions sont donc maintenant incluses "en dur" dans le code, 
# et ne suivront plus les mises à jours éventuelles du package soundecology. A surveiller en cas de bug...
# 5.09 ajout d'un paramètre duree_min (ici et dans calcul_Indices_Acoustiques_batch) pour permettre
# l'analyse de fichiers de petite taille (snipets) sans perdre la possibilité d'ajouter 
# une sécurité contre les fichiers corrompus

library(seewave)
library (tuneR)
library(soundecology)

#----------------------------------------------------------------------------------------------------------
# tdelattre : making bioacoustic function from "soundecology" package SILENT
# version of the function extracted on 28.05.24

bioacSilent <-
function (soundfile, min_freq = 2000, max_freq = 8000, fft_w = 512) 
{
  if (is.numeric(as.numeric(min_freq))) {
    min_freq <- as.numeric(min_freq)
  }
  else {
    stop(" min_freq is not a number.")
  }
  if (is.numeric(as.numeric(max_freq))) {
    max_freq <- as.numeric(max_freq)
  }
  else {
    stop(" max_freq is not a number.")
  }
  if (is.numeric(as.numeric(fft_w))) {
    fft_w <- as.numeric(fft_w)
  }
  else {
    stop(" fft_w is not a number.")
  }
  samplingrate <- soundfile@samp.rate
  freq_per_row = 10
  wlen = samplingrate/freq_per_row
  nyquist_freq <- samplingrate/2
  if (max_freq > nyquist_freq) {
    cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", 
              nyquist_freq, "Hz. But the script was set to measure up to ", 
              max_freq, "Hz.\n\n", sep = ""))
  }
  if (soundfile@stereo == TRUE) {
    #cat("\n This is a stereo file. Results will be given for each channel.\n")  ##tdelattre 28.05.24 making bioac silent
    left <- channel(soundfile, which = c("left"))
    right <- channel(soundfile, which = c("right"))
    #cat("\n Calculating index. Please wait... \n\n")  ##tdelattre 28.05.24 making bioac silent
    spec_left <- spectro(left, f = samplingrate, wl = fft_w, 
                         plot = FALSE, dB = "max0")$amp
    spec_right <- spectro(right, f = samplingrate, wl = fft_w, 
                          plot = FALSE, dB = "max0")$amp
    rm(left, right)
    specA_left <- apply(spec_left, 1, meandB)
    specA_right <- apply(spec_right, 1, meandB)
    rows_width = length(specA_left)/nyquist_freq
    min_row = min_freq * rows_width
    max_row = max_freq * rows_width
    specA_left_segment <- specA_left[min_row:max_row]
    specA_right_segment <- specA_right[min_row:max_row]
    freq_range <- max_freq - min_freq
    freqs <- seq(from = min_freq, to = max_freq, length.out = length(specA_left_segment))
    specA_left_segment_normalized <- specA_left_segment - 
      min(specA_left_segment)
    specA_right_segment_normalized <- specA_right_segment - 
      min(specA_right_segment)
    left_area <- sum(specA_left_segment_normalized * rows_width)
    right_area <- sum(specA_right_segment_normalized * rows_width)
    cat("  Bioacoustic Index:\n")
    cat("   Left channel: ")
    cat(left_area)
    cat("\n   Right channel: ")
    cat(right_area)
    cat("\n\n")
  }
  else {
    #cat("\n This is a mono file.\n")  ##tdelattre 28.05.24 making bioac silent
    left <- channel(soundfile, which = c("left"))
    #cat("\n Calculating index. Please wait... \n\n")  ##tdelattre 28.05.24 making bioac silent
    spec_left <- spectro(left, f = samplingrate, wl = fft_w, 
                         plot = FALSE, dB = "max0")$amp
    rm(left)
    specA_left <- apply(spec_left, 1, meandB)
    rows_width = length(specA_left)/nyquist_freq
    min_row = min_freq * rows_width
    max_row = max_freq * rows_width
    specA_left_segment <- specA_left[min_row:max_row]
    freq_range <- max_freq - min_freq
    freqs <- seq(from = min_freq, to = max_freq, length.out = length(specA_left_segment))
    specA_left_segment_normalized <- specA_left_segment - 
      min(specA_left_segment)
    left_area <- sum(specA_left_segment_normalized * rows_width)
    #cat("  Bioacoustic Index: ")  ##tdelattre 28.05.24 making bioac silent
    #cat(left_area)                ##tdelattre 28.05.24 making bioac silent
    #cat("\n\n")                   ##tdelattre 28.05.24 making bioac silent
    right_area <- NA
  }
  invisible(list(left_area = left_area, right_area = right_area))
}

#-----------------------------------------------------------------------------------------------------------
# tdelattre : making acoustic_complexity function from "soundecology" package SILENT
# version of the function extracted on 28.05.24

acoustic_complexity_Silent<-
function (soundfile, min_freq = NA, max_freq = NA, j = 5, fft_w = 512) 
{
  if (is.na(max_freq)) {
    max_freq <- soundfile@samp.rate/2
    cat(paste("\n max_freq not set, using value of:", max_freq, 
              "\n\n"))
  }
  if (is.na(min_freq)) {
    min_freq <- 0
    cat(paste("\n min_freq not set, using value of:", min_freq, 
              "\n\n"))
  }
  if (is.numeric(as.numeric(min_freq))) {
    min_freq <- as.numeric(min_freq)
  }
  else {
    stop(" min_freq is not a number.")
  }
  if (is.numeric(as.numeric(max_freq))) {
    max_freq <- as.numeric(max_freq)
  }
  else {
    stop(" max_freq is not a number.")
  }
  if (is.numeric(as.numeric(j))) {
    j <- as.numeric(j)
  }
  else {
    stop(" j is not a number.")
  }
  if (is.numeric(as.numeric(fft_w))) {
    fft_w <- as.numeric(fft_w)
  }
  else {
    stop(" fft_w is not a number.")
  }
  get_d <- function(spectrum, freq_row, min_col, max_col) {
    D = 0
    for (k in min_col:(max_col - 1)) {
      D = D + abs(spectrum[freq_row, k] - spectrum[freq_row, 
                                                   k + 1])
    }
    return(D)
  }
  samplingrate <- soundfile@samp.rate
  duration <- length(soundfile@left)/soundfile@samp.rate
  nyquist_freq <- (samplingrate/2)
  if (max_freq > nyquist_freq) {
    cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", 
              nyquist_freq, "Hz. But the script was set to measure up to ", 
              max_freq, "Hz. The value of max_freq was changed to ", 
              nyquist_freq, ".\n\n", sep = ""))
    max_freq <- nyquist_freq
  }
  wlen = fft_w
  if (soundfile@stereo == TRUE) {
    cat("\n This is a stereo file. Results will be given for each channel.\n")
    left <- channel(soundfile, which = c("left"))
    right <- channel(soundfile, which = c("right"))
    cat("\n Calculating index. Please wait... \n\n")
    spec_left <- spectro(left, f = samplingrate, wl = wlen, 
                         plot = FALSE, norm = TRUE, dB = NULL, scale = FALSE, 
                         wn = "hamming")
    specA_left <- spec_left$amp
    min_freq1k = min_freq/1000
    max_freq1k = max_freq/1000
    which_min_freq <- which(abs(spec_left$freq - min_freq1k) == 
                              min(abs(spec_left$freq - min_freq1k)))
    which_max_freq <- which(abs(spec_left$freq - max_freq1k) == 
                              min(abs(spec_left$freq - max_freq1k)))
    if (which_min_freq < 1) {
      which_min_freq = 1
    }
    if (which_max_freq > dim(specA_left)[1]) {
      which_max_freq = dim(specA_left)[1] - 1
    }
    specA_left <- spec_left$amp[which_min_freq:which_max_freq, 
    ]
    rm(spec_left)
    spec_right <- spectro(right, f = samplingrate, wl = wlen, 
                          plot = FALSE, norm = TRUE, dB = NULL, scale = FALSE, 
                          wn = "hamming")
    specA_right <- spec_right$amp[which_min_freq:which_max_freq, 
    ]
    rm(spec_right)
    rm(left, right)
    specA_rows <- dim(specA_left)[1]
    specA_cols <- dim(specA_left)[2]
    fl <- rep(NA, specA_rows)
    delta_fl <- (max_freq - min_freq)/specA_rows
    delta_tk <- (length(soundfile@left)/soundfile@samp.rate)/specA_cols
    no_j <- floor(duration/j)
    I_per_j <- floor(j/delta_tk)
    ACI_left_vals <- rep(NA, no_j)
    ACI_fl_left_vector <- rep(NA, no_j)
    ACI_left_matrix <- data.frame(matrix(NA, nrow = specA_rows, 
                                         ncol = no_j))
    ACI_right_vals <- rep(NA, no_j)
    ACI_fl_right_vector <- rep(NA, no_j)
    ACI_right_matrix <- data.frame(matrix(NA, nrow = specA_rows, 
                                          ncol = no_j))
    for (q_index in 1:specA_rows) {
      for (j_index in 1:no_j) {
        min_col <- j_index * I_per_j - I_per_j + 1
        max_col <- j_index * I_per_j
        D <- get_d(specA_left, q_index, min_col, max_col)
        sum_I <- sum(specA_left[q_index, min_col:max_col])
        ACI_left_vals[j_index] <- D/sum_I
        ACI_left_matrix[q_index, j_index] <- D/sum_I
      }
      ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
    }
    ACI_tot_left <- sum(ACI_fl_left_vector)
    for (q_index in 1:specA_rows) {
      for (j_index in 1:no_j) {
        min_col <- j_index * I_per_j - I_per_j + 1
        max_col <- j_index * I_per_j
        D <- get_d(specA_right, q_index, min_col, max_col)
        sum_I <- sum(specA_right[q_index, min_col:max_col])
        ACI_right_vals[j_index] <- D/sum_I
        ACI_right_matrix[q_index, j_index] <- D/sum_I
      }
      ACI_fl_right_vector[q_index] <- sum(ACI_right_vals)
    }
    ACI_tot_right <- sum(ACI_fl_right_vector)
    ACI_tot_left_by_min <- round((ACI_tot_left/duration) * 
                                   60, 2)
    ACI_tot_right_by_min <- round((ACI_tot_right/duration) * 
                                    60, 2)
    cat(paste("  Acoustic Complexity Index (total):\n", "   Left channel: ", 
              sep = ""))
    cat(ACI_tot_left)
    cat(paste("\n", "   Right channel: ", sep = ""))
    cat(ACI_tot_right)
    cat("\n\n")
    if (duration > 60) {
      cat(paste("  Acoustic Complexity Index (by minute):\n", 
                "   Left channel: ", sep = ""))
      cat(ACI_tot_left_by_min)
      cat(paste("\n", "   Right channel: ", sep = ""))
      cat(ACI_tot_right_by_min)
      cat("\n\n")
    }
  }
  else {
    #cat("\n This is a mono file.\n")
    left <- channel(soundfile, which = c("left"))
    #cat("\n Calculating index. Please wait... \n\n")
    spec_left <- spectro(left, f = samplingrate, wl = wlen, 
                         plot = FALSE, norm = TRUE, dB = NULL, scale = FALSE, 
                         wn = "hamming")
    specA_left <- spec_left$amp
    min_freq1k = min_freq/1000
    max_freq1k = max_freq/1000
    which_min_freq <- which(abs(spec_left$freq - min_freq1k) == 
                              min(abs(spec_left$freq - min_freq1k)))
    which_max_freq <- which(abs(spec_left$freq - max_freq1k) == 
                              min(abs(spec_left$freq - max_freq1k)))
    specA_left <- specA_left[which_min_freq:which_max_freq, 
    ]
    rm(spec_left)
    rm(left)
    specA_rows <- dim(specA_left)[1]
    specA_cols <- dim(specA_left)[2]
    fl <- rep(NA, specA_rows)
    delta_fl <- (max_freq - min_freq)/specA_rows
    delta_tk <- (length(soundfile@left)/soundfile@samp.rate)/specA_cols
    no_j <- floor(duration/j)
    I_per_j <- floor(j/delta_tk)
    ACI_left_vals <- rep(NA, no_j)
    ACI_fl_left_vector <- rep(NA, no_j)
    ACI_left_matrix <- data.frame(matrix(NA, nrow = specA_rows, 
                                         ncol = no_j))
    ACI_right_vals <- rep(NA, no_j)
    ACI_fl_right_vector <- rep(NA, no_j)
    ACI_right_matrix <- data.frame(matrix(NA, nrow = specA_rows, 
                                          ncol = no_j))
    for (q_index in 1:specA_rows) {
      for (j_index in 1:no_j) {
        min_col <- j_index * I_per_j - I_per_j + 1
        max_col <- j_index * I_per_j
        D <- get_d(specA_left, q_index, min_col, max_col)
        sum_I <- sum(specA_left[q_index, min_col:max_col])
        ACI_left_vals[j_index] <- D/sum_I
        ACI_left_matrix[q_index, j_index] <- D/sum_I
      }
      ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
    }
    ACI_tot_left <- sum(ACI_fl_left_vector)
    ACI_tot_left_by_min <- round((ACI_tot_left/duration) * 
                                   60, 2)
    ACI_tot_right <- NA
    ACI_tot_right_by_min <- NA
    #cat("  Acoustic Complexity Index (total): ")
    #cat(ACI_tot_left)
    #cat("\n\n")
    if (duration > 60) {
      cat("  Acoustic Complexity Index (by minute): ")
      cat(ACI_tot_left_by_min)
      cat("\n\n")
    }
  }
  invisible(list(AciTotAll_left = ACI_tot_left, AciTotAll_right = ACI_tot_right, 
                 AciTotAll_left_bymin = ACI_tot_left_by_min, AciTotAll_right_bymin = ACI_tot_right_by_min, 
                 aci_fl_left_vals = ACI_fl_left_vector, aci_fl_right_vals = ACI_fl_right_vector, 
                 aci_left_matrix = ACI_left_matrix, aci_right_matrix = ACI_right_matrix))
}

# tdelattre : making XX function from "soundecology" package SILENT
# version of the function extracted on 28.05.24
#-----------------------------------------------------------------------------------------------------------
# tdelattre : making ndsi function from "soundecology" package SILENT
# version of the function extracted on 28.05.24
#-------------------------------------------------------------------------------------------------------------------------------------
ndsi_Silent <-
function (soundfile, fft_w = 1024, anthro_min = 1000, anthro_max = 2000, 
          bio_min = 2000, bio_max = 11000) 
{
  if (is.numeric(as.numeric(fft_w))) {
    fft_w <- as.numeric(fft_w)
  }
  else {
    stop(" fft_w is not a number.")
  }
  if (is.numeric(as.numeric(anthro_min))) {
    anthro_min <- as.numeric(anthro_min)
  }
  else {
    stop(" anthro_min is not a number.")
  }
  if (is.numeric(as.numeric(anthro_max))) {
    anthro_max <- as.numeric(anthro_max)
  }
  else {
    stop(" anthro_max is not a number.")
  }
  if (is.numeric(as.numeric(bio_min))) {
    bio_min <- as.numeric(bio_min)
  }
  else {
    stop(" bio_min is not a number.")
  }
  if (is.numeric(as.numeric(bio_max))) {
    bio_max <- as.numeric(bio_max)
  }
  else {
    stop(" bio_max is not a number.")
  }
  hz_interval = anthro_max - anthro_min
  samplingrate <- soundfile@samp.rate
  duration <- length(soundfile@left)/soundfile@samp.rate
  nyquist_freq <- (samplingrate/2)
  if (bio_max > nyquist_freq) {
    stop(paste("The maximum frequency of biophony (", bio_max, 
               " Hz) can not be higher than the maximum frequency of the file (", 
               nyquist_freq, " Hz)\n\n Change the value of bio_max to less than ", 
               nyquist_freq, "\n\n", sep = ""))
  }
  if (anthro_max > bio_min) {
    stop(paste("The maximum frequency of anthrophony (", 
               anthro_max, " Hz) can not be higher than the minimum frequency of biophony (", 
               bio_min, " Hz)\n\n Change the value of anthro_max to equal or less than bio_min\n\n", 
               sep = ""))
  }
  if (anthro_max < anthro_min) {
    stop(paste("The minimum frequency of anthrophony (", 
               anthro_min, " Hz) can not be higher than the maximum frequency of anthrophony (", 
               anthro_max, " Hz)\n\n Change the value of anthro_min to less than anthro_max\n\n", 
               sep = ""))
  }
  if (bio_max < bio_min) {
    stop(paste("The minimum frequency of biophony (", bio_min, 
               " Hz) can not be higher than the maximum frequency of biophony (", 
               bio_max, " Hz)\n\n Change the value of anthro_min to less than anthro_max\n\n", 
               sep = ""))
  }
  if (soundfile@stereo == TRUE) {
    cat("\n This is a stereo file. Results will be given for each channel.\n")
    left <- channel(soundfile, which = c("left"))
    right <- channel(soundfile, which = c("right"))
    rm(soundfile)
    cat("\n Calculating index. Please wait... \n")
    left1 <- as.vector(cutw(left, from = 0, to = length(left@left)/left@samp.rate))
    left2 <- data.frame(matrix(NA, nrow = samplingrate, ncol = floor(duration)))
    for (i in 0:(floor(duration) - 1)) {
      j <- i + 1
      start1 <- (i * samplingrate) + 1
      end <- start1 + samplingrate - 1
      left2[, j] <- left1[start1:end]
    }
    left3 <- data.frame(matrix(NA, nrow = fft_w/2, ncol = floor(duration)))
    left4 <- apply(left2, 2, oce::pwelch, fs = samplingrate, nfft = fft_w, 
                   plot = FALSE)
    for (i in 1:floor(duration)) {
      left3[, i] <- left4[[i]]$spec
    }
    specA_left <- apply(left3, 1, mean)
    specA_rows <- length(specA_left)
    freq_per_row <- specA_rows/nyquist_freq
    anthro_vals_range <- anthro_max - anthro_min
    bio_vals_range <- bio_max - bio_min
    bio_bins <- round(bio_vals_range/hz_interval)
    anthro_bins <- rep(NA, round(anthro_vals_range/hz_interval))
    bio_bins <- rep(NA, round(bio_vals_range/hz_interval))
    anthro_min_row <- round(anthro_min * freq_per_row)
    anthro_max_row <- round(anthro_max * freq_per_row)
    bio_step_range <- freq_per_row * (bio_vals_range/length(bio_bins))
    bio_min_row <- round(bio_min * freq_per_row)
    bio_max_row <- bio_min_row + bio_step_range
    for (i in 1:length(anthro_bins)) {
      anthro_bins[i] <- pracma::trapz(specA_left[anthro_min_row:anthro_max_row])
    }
    for (i in 1:length(bio_bins)) {
      if (bio_max_row >= specA_rows) {
        bio_max_row <- specA_rows
      }
      bio_bins[i] <- pracma::trapz(specA_left[bio_min_row:bio_max_row])
      bio_min_row <- bio_min_row + bio_step_range
      bio_max_row <- bio_max_row + bio_step_range
    }
    freqbins <- rep(NA, sum(length(anthro_bins), length(bio_bins)))
    freqbins <- c(anthro_bins, bio_bins)
    freqbins = freqbins/norm(as.matrix(freqbins), "F")
    freqbins.SumAll <- sum(freqbins)
    freqbins.SumBio <- sum(freqbins[2:length(freqbins)])
    freqbins.Anthro <- freqbins[1]
    NDSI_left <- (freqbins.SumBio - freqbins.Anthro)/(freqbins.SumBio + 
                                                        freqbins.Anthro)
    biophony_left <- freqbins.SumBio
    anthrophony_left <- freqbins.Anthro
    right1 <- as.vector(cutw(right, from = 0, to = length(right@left)/right@samp.rate))
    right2 <- data.frame(matrix(NA, nrow = samplingrate, 
                                ncol = floor(duration)))
    for (i in 0:(floor(duration) - 1)) {
      j <- i + 1
      start1 <- (i * samplingrate) + 1
      end <- start1 + samplingrate - 1
      right2[, j] <- right1[start1:end]
    }
    right3 <- data.frame(matrix(NA, nrow = fft_w/2, ncol = floor(duration)))
    right4 <- apply(right2, 2,oce::pwelch, fs = samplingrate, 
                    nfft = fft_w, plot = FALSE)
    for (i in 1:floor(duration)) {
      right3[, i] <- right4[[i]]$spec
    }
    specA_right <- apply(right3, 1, mean)
    specA_rows <- length(specA_right)
    freq_per_row <- specA_rows/nyquist_freq
    anthro_vals_range <- anthro_max - anthro_min
    bio_vals_range <- bio_max - bio_min
    bio_bins <- round(bio_vals_range/hz_interval)
    anthro_bins <- rep(NA, round(anthro_vals_range/hz_interval))
    bio_bins <- rep(NA, round(bio_vals_range/hz_interval))
    anthro_min_row <- round(anthro_min * freq_per_row)
    anthro_max_row <- round(anthro_max * freq_per_row)
    bio_step_range <- freq_per_row * (bio_vals_range/length(bio_bins))
    bio_min_row <- round(bio_min * freq_per_row)
    bio_max_row <- bio_min_row + bio_step_range
    for (i in 1:length(anthro_bins)) {
      anthro_bins[i] <- pracma::trapz(specA_right[anthro_min_row:anthro_max_row])
    }
    for (i in 1:length(bio_bins)) {
      if (bio_max_row >= specA_rows) {
        bio_max_row <- specA_rows
      }
      bio_bins[i] <- pracma::trapz(specA_right[bio_min_row:bio_max_row])
      bio_min_row <- bio_min_row + bio_step_range
      bio_max_row <- bio_max_row + bio_step_range
    }
    freqbins <- rep(NA, sum(length(anthro_bins), length(bio_bins)))
    freqbins <- c(anthro_bins, bio_bins)
    freqbins = freqbins/norm(as.matrix(freqbins), "F")
    freqbins.SumAll <- sum(freqbins)
    freqbins.SumBio <- sum(freqbins[2:length(freqbins)])
    freqbins.Anthro <- freqbins[1]
    NDSI_right <- (freqbins.SumBio - freqbins.Anthro)/(freqbins.SumBio + 
                                                         freqbins.Anthro)
    biophony_right <- freqbins.SumBio
    anthrophony_right <- freqbins.Anthro
    cat("  Normalized Difference Soundscape Index:\n")
    cat("\n   Left channel: ")
    cat(NDSI_left)
    cat("\n")
    cat("   Right channel: ")
    cat(NDSI_right)
    cat("\n\n")
  }
  else {
    #cat("\n This is a mono file.\n")
    left <- channel(soundfile, which = c("left"))
    rm(soundfile)
    #cat("\n Calculating index. Please wait... \n\n")
    left1 <- as.vector(cutw(left, from = 0, to = length(left@left)/left@samp.rate))
    left2 <- data.frame(matrix(NA, nrow = samplingrate, ncol = floor(duration)))
    for (i in 0:(floor(duration) - 1)) {
      j <- i + 1
      start1 <- (i * samplingrate) + 1
      end <- start1 + samplingrate - 1
      left2[, j] <- left1[start1:end]
    }
    left3 <- data.frame(matrix(NA, nrow = fft_w/2, ncol = floor(duration)))
    left4 <- apply(left2, 2,oce::pwelch, fs = samplingrate, nfft = fft_w, 
                   plot = FALSE)
    for (i in 1:floor(duration)) {
      left3[, i] <- left4[[i]]$spec
    }
    specA_left <- apply(left3, 1, mean)
    specA_rows <- length(specA_left)
    freq_per_row <- specA_rows/nyquist_freq
    anthro_vals_range <- anthro_max - anthro_min
    bio_vals_range <- bio_max - bio_min
    bio_bins <- round(bio_vals_range/hz_interval)
    anthro_bins <- rep(NA, round(anthro_vals_range/hz_interval))
    bio_bins <- rep(NA, round(bio_vals_range/hz_interval))
    anthro_min_row <- round(anthro_min * freq_per_row)
    anthro_max_row <- round(anthro_max * freq_per_row)
    bio_step_range <- freq_per_row * (bio_vals_range/length(bio_bins))
    bio_min_row <- round(bio_min * freq_per_row)
    bio_max_row <- bio_min_row + bio_step_range
    for (i in 1:length(anthro_bins)) {
      anthro_bins[i] <- pracma::trapz(specA_left[anthro_min_row:anthro_max_row])
    }
    for (i in 1:length(bio_bins)) {
      if (bio_max_row >= specA_rows) {
        bio_max_row <- specA_rows
      }
      bio_bins[i] <- pracma::trapz(specA_left[bio_min_row:bio_max_row])
      bio_min_row <- bio_min_row + bio_step_range
      bio_max_row <- bio_max_row + bio_step_range
    }
    freqbins <- rep(NA, sum(length(anthro_bins), length(bio_bins)))
    freqbins <- c(anthro_bins, bio_bins)
    freqbins = freqbins/norm(as.matrix(freqbins), "F")
    freqbins.SumAll <- sum(freqbins)
    freqbins.SumBio <- sum(freqbins[2:length(freqbins)])
    freqbins.Anthro <- freqbins[1]
    NDSI_left <- (freqbins.SumBio - freqbins.Anthro)/(freqbins.SumBio + 
                                                        freqbins.Anthro)
    biophony_left <- freqbins.SumBio
    anthrophony_left <- freqbins.Anthro
    biophony_right <- NA
    anthrophony_right <- NA
    NDSI_right = NA
    #cat("  Normalized Difference Soundscape Index: ")
    #cat(NDSI_left)
    #cat("\n\n")
  }
  invisible(list(ndsi_left = NDSI_left, ndsi_right = NDSI_right, 
                 biophony_left = biophony_left, anthrophony_left = anthrophony_left, 
                 biophony_right = biophony_right, anthrophony_right = anthrophony_right))
}
# tdelattre : making XX function from "soundecology" package SILENT
# version of the function extracted on 28.05.24
#-------------------------------------------------------------------------------------------------------------------------------------
# tdelattre : making acoustic_diversity function from "soundecology" package SILENT
# version of the function extracted on 28.05.24
#-------------------------------------------------------------------------------------------------------------------------------------
acoustic_diversity_Silent<-
function (soundfile, max_freq = 10000, db_threshold = -50, freq_step = 1000, 
          shannon = TRUE) 
{
  db_threshold <- as.numeric(db_threshold)
  if (is.numeric(as.numeric(max_freq))) {
    max_freq <- as.numeric(max_freq)
  }
  else {
    stop(" max_freq is not a number.")
  }
  if (is.numeric(as.numeric(db_threshold))) {
    db_threshold <- as.numeric(db_threshold)
  }
  else {
    stop(" db_threshold is not a number.")
  }
  if (is.numeric(as.numeric(freq_step))) {
    freq_step <- as.numeric(freq_step)
  }
  else {
    stop(" freq_step is not a number.")
  }
  getscore <- function(spectrum, minf, maxf, db, freq_row) {
    miny <- round((minf)/freq_row)
    maxy <- round((maxf)/freq_row)
    subA = spectrum[miny:maxy, ]
    index1 <- length(subA[subA > db])/length(subA)
    return(index1)
  }
  samplingrate <- soundfile@samp.rate
  nyquist_freq <- samplingrate/2
  freq_per_row = 10
  wlen = samplingrate/freq_per_row
  if (wlen%%2 == 1) {
    wlen <- wlen + 1
  }
  if (soundfile@stereo == TRUE) {
    cat("\n This is a stereo file. Results will be given for each channel.\n")
    left <- channel(soundfile, which = c("left"))
    right <- channel(soundfile, which = c("right"))
    rm(soundfile)
    cat("\n Calculating index. Please wait... \n\n")
    specA_left <- spectro(left, f = samplingrate, wl = wlen, 
                          plot = FALSE)$amp
    specA_right <- spectro(right, f = samplingrate, wl = wlen, 
                           plot = FALSE)$amp
    rm(left, right)
    if (max_freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", 
                nyquist_freq, "Hz. But the script was set to measure up to ", 
                max_freq, "Hz. The value of max_freq was changed to ", 
                nyquist_freq, ".\n\n", sep = ""))
      max_freq <- nyquist_freq
    }
    Freq <- seq(from = 0, to = max_freq - freq_step, by = freq_step)
    Score <- rep(NA, length(Freq))
    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_left, Freq[j], (Freq[j] + 
                                                  freq_step), db_threshold, freq_per_row)
    }
    left_vals = Score
    Score1 = 0
    for (i in 1:length(Freq)) {
      Score1 = Score1 + (Score[i] * log(Score[i] + 1e-07))
    }
    Score_left = (-(Score1))/length(Freq)
    Shannon_left <- vegan::diversity(Score, index = "shannon")
    Score <- rep(NA, length(Freq))
    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_right, Freq[j], (Freq[j] + 
                                                   freq_step), db_threshold, freq_per_row)
    }
    right_vals = Score
    Score1 = 0
    for (i in 1:length(Freq)) {
      Score1 = Score1 + (Score[i] * log(Score[i] + 1e-07))
    }
    Score_right = (-(Score1))/length(Freq)
    Shannon_right <- vegan::diversity(Score, index = "shannon")
    left_bandvals_return <- rep(NA, length(Freq))
    right_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    right_bandrange_return <- rep(NA, length(Freq))
    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = round(left_vals[j], 6)
      right_bandvals_return[j] = round(right_vals[j], 6)
      left_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + 
                                                        freq_step), " Hz", sep = "")
      right_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + 
                                                         freq_step), " Hz", sep = "")
    }
    for (j in seq(length(Freq), 1, by = -1)) {
      this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), 
                             "", sep = "")
      this_row_size <- nchar(this_row_name)
      this_row_space <- 17 - this_row_size
      this_row_spaces = ""
      for (f in seq(1, this_row_space, by = 1)) {
        this_row_spaces = paste(this_row_spaces, " ", 
                                sep = "")
      }
    }
    for (j in seq(length(Freq), 1, by = -1)) {
      this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), 
                             "", sep = "")
      this_row_size <- nchar(this_row_name)
      this_row_space <- 17 - this_row_size
      this_row_spaces = ""
      for (f in seq(1, this_row_space, by = 1)) {
        this_row_spaces = paste(this_row_spaces, " ", 
                                sep = "")
      }
    }
    if (shannon == TRUE) {
      left_adi_return = round(Shannon_left, 6)
      right_adi_return = round(Shannon_right, 6)
    }
    else {
      left_adi_return = round(Score_left, 6)
      right_adi_return = round(Score_right, 6)
    }
    cat("  Acoustic Diversity Index: \n")
    cat(paste("   Left channel: ", left_adi_return, "\n", 
              sep = ""))
    cat(paste("   Right channel: ", right_adi_return, "\n", 
              sep = ""))
  }
  else {
    #cat("\n This is a mono file.\n")
    #cat("\n Calculating index. Please wait... \n\n")
    specA_left <- spectro(soundfile, f = samplingrate, wl = wlen, 
                          plot = FALSE)$amp
    rm(soundfile)
    if (max_freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", 
                nyquist_freq, "Hz. But the script was set to measure up to ", 
                max_freq, "Hz. The value of max_freq was changed to ", 
                nyquist_freq, ".\n\n", sep = ""))
      max_freq <- nyquist_freq
    }
    Freq <- seq(from = 0, to = max_freq - freq_step, by = freq_step)
    Score <- rep(NA, length(Freq))
    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_left, Freq[j], (Freq[j] + 
                                                  freq_step), db_threshold, freq_per_row)
    }
    left_vals = Score
    Score1 = 0
    for (i in 1:length(Freq)) {
      Score1 = Score1 + (Score[i] * log(Score[i] + 1e-07))
    }
    Score_left = (-(Score1))/length(Freq)
    Shannon_left <- vegan::diversity(Score, index = "shannon")
    Shannon_right <- NA
    left_bandvals_return <- rep(NA, length(Freq))
    right_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    right_bandrange_return <- rep(NA, length(Freq))
    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = round(left_vals[j], 6)
      left_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + 
                                                        freq_step), " Hz", sep = "")
    }
    for (j in seq(length(Freq), 1, by = -1)) {
      this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), 
                             "", sep = "")
      this_row_size <- nchar(this_row_name)
      this_row_space <- 17 - this_row_size
      this_row_spaces = ""
      for (f in seq(1, this_row_space, by = 1)) {
        this_row_spaces = paste(this_row_spaces, " ", 
                                sep = "")
      }
    }
    #cat("  Acoustic Diversity Index: ")
    right_adi_return = NA
    if (shannon == TRUE) {
      #cat(paste(round(Shannon_left, 6), "\n", sep = ""))
      left_adi_return = round(Shannon_left, 6)
    }
    else {
      #cat(paste(round(Score_left, 6), "\n", sep = ""))
      left_adi_return = round(Score_left, 6)
    }
  }
  invisible(list(adi_left = left_adi_return, adi_right = right_adi_return, 
                 left_band_values = left_bandvals_return, right_band_values = right_bandvals_return, 
                 left_bandrange_values = left_bandrange_return, right_bandrange_values = right_bandrange_return))
}
# tdelattre : making XX function from "soundecology" package SILENT
# version of the function extracted on 28.05.24
#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# slight modification of the fpeak function from the "seewave" package.
#-------------------------------------------------------------------------------------------------------------------------------------
fpeaksFlat <-
  function (spec,
            f = NULL,
            nmax = NULL,
            amp = NULL,
            freq = NULL,
            threshold = NULL,
            plot = TRUE,
            title = TRUE,
            xlab = "Frequency (kHz)",
            ylab = "Amplitude",
            labels = TRUE,
            legend = TRUE,
            collab = "red",
            ...)
  {
    if (is.matrix(spec)) {
      if (ncol(spec) != 2)
        stop(
          "If 'spec' is a numeric matrix it should be a two column matrix with the first colum describing the frequency x-axis and the second column describing the amplitude y-axis"
        )
      N <- nrow(spec)
    }
    if (is.vector(spec)) {
      N <- length(spec)
      if (is.null(f)) {
        stop(
          "If 'spec' is a numeric vector describing the amplitude only, the sampling frequency 'f' of the original signal should be provided (for instance (f = 44100)"
        )
      }
      if (!is.null(f) && !is.na(f)) {
        spec <- cbind(seq(f / (N * 2), f / 2, length = N) / 1000,
                      spec)
      }
      if (!is.null(f) && is.na(f)) {
        spec <- cbind(1:N, spec)
        plot <- FALSE
      }
    }
    flat <- round(N / 20)
    spec.tmp <- c(spec[, 2], rep(NA, flat))
    for (i in 1:(N - (flat + 1))) {
      ref <- spec.tmp[i]
      for (j in 1:flat) {
        if (spec.tmp[i + j] == ref) {
          spec.tmp[i + j] <- spec.tmp[i + j] + 1e-05 *
            spec.tmp[i + j]
        }
      }
    }
    spec <- cbind(spec[, 1], spec.tmp[1:N])
    sym <- discrets(spec[, 2], symb = 5, collapse = FALSE)
    
    # si ya des flats
    specF <- spec
    for (mod in 1:length(specF[, 2]))
    {
      if (specF[mod, 2] == 0) {
        specF[mod, 2] <- specF[mod, 2] + 0.000001 * mod
      }
    }
    sym2 <- discrets(specF[, 2], symb = 5, collapse = FALSE)
    #
    
    if (sym2[1] == "I")
      sym2[1] <- "T"
    if (sym2[1] == "P")
      sym2[1] <- "D"
    
    sym2 <- c(NA, sym2, NA)
    peaks <- which(sym2 == "P")
    valleys <- which(sym2 == "T")
    
    n <- length(peaks)
    if (n == 0) {
      res <- NA
      plot <- FALSE
    }    else {
      if (!is.null(amp) | !is.null(nmax)) {
        diffvp <- diffpv <- numeric(n)
        for (i in 1:n) {
          v <- specF[valleys[i], 2]
          p <- specF[peaks[i], 2]
          vv <- specF[valleys[i + 1], 2]
          diffvp[i] <- p - v
          diffpv[i] <- p - vv
        }
      }
      if (!is.null(nmax) && n != 0) {
        if (!is.null(amp) | !is.null(freq) | !is.null(threshold)) {
          cat(
            "Caution! The argument 'nmax' overrides the arguments 'amp', 'freq', and 'threshold'"
          )
        }
        if (n < nmax) {
          cat(paste("There are", n, "peaks only (< nmax ="),
              nmax, ")")
        }
        if (nmax == 1) {
          tmp <- specF[peaks, , drop = FALSE]
          res <- tmp[which.max(tmp[, 2]), , drop = FALSE]
        }            else {
          alt <- cbind(peaks, diffvp, diffpv)
          leftorder <- alt[order(-alt[, 2]), , drop = FALSE]
          rightorder <- alt[order(-alt[, 3]), , drop = FALSE]
          left <- leftorder[, 1]
          right <- rightorder[, 1]
          l <- 0
          i <- 1
          while (l[i] < nmax) {
            comp <- left[1:i] %in% right[1:i]
            l <- c(l, length(comp[comp == TRUE]))
            i <- i + 1
          }
          peaks0 <- left[1:(i - 1)]
          if (l[i] > nmax) {
            error <- l[i] - nmax
            peaks0 <- peaks0[1:(length(peaks0) - error)]
          }
          peaks <- peaks0[comp]
          res <- matrix(na.omit(specF[peaks,]), nc = 2)
          colnames(res) <- c("freq", "amp")
        }
      } else {
        if (!is.null(amp)) {
          if (length(amp) != 2)
            stop("The length of 'amp' should equal to 2.")
          for (i in 1:n) {
            if (!is.na(diffvp[i]) && !is.na(diffpv[i]) &&
                diffvp[i] > 0 && diffpv[i] > 0 && diffvp[i] >=
                amp[1] && diffpv[i] >= amp[2])
              peaks[i] <- peaks[i]
            else
              peaks[i] <- NA
          }
        }
        if (!is.null(freq)) {
          freq <- freq / 1000
          diffpeak <- numeric(n - 1)
          for (i in 1:(n - 1)) {
            peak1 <- specF[peaks[i], 1]
            peak2 <- specF[peaks[i + 1], 1]
            diffpeak[i] <- peak2 - peak1
            if (!is.na(diffpeak[i]) && diffpeak[i] <= freq)
              if (specF[peaks[i + 1], 2] > specF[peaks[i],
                                                 2])
              {
                peaks[i] <- NA
              }  else
                peaks[i + 1] <- NA
          }
        }
        if (!is.null(f) && is.na(f)) {
          res <- peaks
        }            else {
          res <- matrix(na.omit(specF[peaks,]), nc = 2)
          colnames(res) <- c("freq", "amp")
        }
        if (!is.null(threshold)) {
          res <- res[res[, 2] > threshold, , drop = FALSE]
        }
      }
    }
    if (plot) {
      plot(
        spec,
        type = "l",
        xlab = xlab,
        ylab = ylab,
        xaxs = "i",
        yaxt = "n",
        ...
      )
      if (title) {
        if (nrow(res) == 1) {
          text.title <- "peak detected"
        }            else {
          text.title <- "peaks detected"
        }
        title(main = paste(nrow(res), text.title))
      }
      points(res, col = collab)
      if (labels & nrow(res) != 0)
        text(res,
             labels = round(res[, 1], 2),
             pos = 3,
             col = collab)
      if (!is.null(threshold)) {
        abline(h = threshold, col = collab, lty = 2)
        mtext(
          paste(threshold),
          side = 2,
          line = 0.5,
          at = threshold,
          las = 1,
          col = collab
        )
      }
      if (legend) {
        if (!is.null(nmax)) {
          text.legend <- paste("nmax=", nmax, sep = "")
        }            else {
          if (is.null(amp)) {
            amp[1] <- amp[2] <- "-"
          }                else
            amp <- round(amp, 2)
          if (is.null(freq)) {
            freq <- "-"
          }
          if (is.null(threshold)) {
            threshold <- "-"
          }
          text.legend <- c(
            paste("amp=", amp[1], "/", amp[2],
                  sep = ""),
            paste("freq=", freq, sep = ""),
            paste("threshold=", threshold, sep = "")
          )
        }
        legend(
          "topright",
          pch = NA,
          legend = text.legend,
          bty = "n",
          text.col = "darkgrey"
        )
      }
      invisible(res)
    }    else
      return(res)
  }



#############################
#############################

AcouIndexAlpha <-
  function(wave,
           mydir,
           stereo = FALSE,
           min_freq = 2000,
           max_freq = 22000,
           min_duration = 57,
           anthro_min = 1000,
           anthro_max = 2000,
           bio_min = 2000,
           bio_max = 12000,
           wl = 512,
           j = 5,
           Bioac = TRUE,
           ACI = TRUE,
           NDSI = TRUE,
           ADI = TRUE,
           NP = TRUE)
  {
    library(seewave)
    library(tuneR)
    library(soundecology)
    
    #ajout tdelattre  pour automatisation sur nos arborescences
    print(paste(mydir, wave, sep = ""))
    wave = readWave(paste(mydir, wave, sep = ""))
    if (seewave::duration(wave) > min_duration) { #pour éviter les fichiers corrompus
      #arguments fpeak
      amp = c(1 / 90, 1 / 90)
      freq = 200
      plotpic = FALSE
      
      f <- wave@samp.rate
      nyquist_freq <- f / 2
      
      left <- channel(wave, which = c("left"))
      spec_left <-
        spectro(
          left,
          f = f,
          wl = wl,
          plot = FALSE,
          dB = "max0"
        )$amp
      specA_left <- apply(spec_left, 1, meandB)
      rows_width = length(specA_left) / nyquist_freq
      min_row = round(min_freq * rows_width)
      max_row = round(max_freq * rows_width)
      specA_left_segment <- specA_left[min_row:max_row]
      freqs <-
        seq(
          from = min_freq,
          to = max_freq,
          length.out = length(specA_left_segment)
        )
      specA_left_segment_normalized <-
        (specA_left_segment - min(specA_left_segment)) / max(specA_left_segment -
                                                               min(specA_left_segment))
      spec_L <- cbind(freqs, specA_left_segment_normalized)
      
      if (wave@stereo == TRUE)
      {
        right <- channel(wave, which = c("right"))
        spec_right <-
          spectro(
            right,
            f = f,
            wl = wl,
            plot = FALSE,
            dB = "max0"
          )$amp
        specA_right <- apply(spec_left, 1, meandB)
        specA_right_segment <- specA_right[min_row:max_row]
        specA_right_segment_normalized <-
          (specA_right_segment - min(specA_right_segment)) / max(specA_right_segment -
                                                                   min(specA_right_segment))
        spec_R <- cbind(freqs, specA_right_segment_normalized)
      }
      
      Table_left <- NULL
      Table_right <- NULL
      
      
      
      ###################################################################
      # Relative avian abondance - (Boelman et al. 2007)
      #calculation of the dB en fonction to the frequency=meanspectrum db using
      #then calculation of the area under the spectrum dB
      if (Bioac == TRUE)
      {
        bioacouMeasure <-
          bioacSilent(wave,
                            min_freq = min_freq,
                            max_freq = max_freq,
                            fft_w = wl)
        Bioac_left <- bioacouMeasure$left_area
        Table_left <- cbind(Table_left, Bioac_left)
        if (wave@stereo == TRUE)
        {
          Bioac_right <- bioacouMeasure$right_area
          Table_right <- cbind(Table_right, Bioac_right)
        }
      }
      
      
      
      
      ###################################################################
      # ACI (Pieretti et al. 2011) ACI() {seewave} A REVOIR A PAROPOS DES FILTRES
      # Issue with the j with short sound# maybe by default =1 ???
      if (ACI == TRUE)
      {
        ACI_measure <-
          acoustic_complexity_Silent(
            wave,
            min_freq = min_freq,
            max_freq = max_freq,
            j = j,
            fft_w = wl
          )
        ACI_left <- ACI_measure$AciTotAll_left
        Table_left <- cbind(Table_left, ACI_left)
        if (wave@stereo == TRUE)
        {
          ACI_right <- ACI_measure$AciTotAll_right
          Table_right <- cbind(Table_right, ACI_right)
        }
      }
      
      
      
      
      
      ###################################################################
      # Normalised difference soundscape index - NDSI (Kasten et al. 2012)
      if (NDSI == TRUE)
      {
        NDSI_measure <-
          ndsi_Silent(
            wave,
            fft_w = wl,
            anthro_min = anthro_min,
            anthro_max = anthro_max,
            bio_min = bio_min,
            bio_max = bio_max
          )
        NDSI_left <- NDSI_measure$ndsi_left
        #NDSI_seewave<-NDSI(soundscapespec(wave))
        Table_left <- cbind(Table_left, NDSI_left)
        
        if (wave@stereo == TRUE)
        {
          NDSI_right <- NDSI_measure$ndsi_right
          Table_right <- cbind(Table_right, NDSI_right)
        }
      }
      
      
      
      
      
      
      
      ###################################################################
      # Acoustic diversity index - ADI (=H') (Pekin et al. 2013)// H' (Villanueva-Riviera et al. 2011)
      if (ADI == TRUE)
      {
        ADI_measure <-
          acoustic_diversity_Silent(
            wave,
            max_freq = max_freq,
            db_threshold = -50,
            freq_step = 1000,
            shannon = TRUE
          )
        ADI_left <- ADI_measure$adi_left
        Table_left <- cbind(Table_left, ADI_left)
        
        if (wave@stereo == TRUE)
        {
          ADI_right <- ADI_measure$adi_right
          Table_right <- cbind(Table_right, ADI_right)
        }
      }
      
      
      
      
      ###################################################################
      # Number of peaks - NP (Gasc et al. 2013b)
      if (NP == TRUE)
      {
        res1_left <- fpeaksFlat(spec_L, plot = F, f)
        pictot_left <- nrow(res1_left)
        
        if (is.null(pictot_left) == FALSE)
        {
          if (pictot_left != 1)
          {
            res2_left <- fpeaksFlat(spec_L,
                                    amp = amp,
                                    freq = freq,
                                    plot = plotpic,
                                    f)
            npic_left <- nrow(res2_left)
          } else{
            npic_left <- 1
          }
        } else{
          npic_left <- 0
        }
        Table_left <- cbind(Table_left, npic_left)
        
        if (wave@stereo == TRUE)
        {
          res1_right <- fpeaksFlat(spec_R, plot = F, f)
          pictot_right <- nrow(res1_right)
          
          if (is.null(pictot_right) == FALSE)
          {
            if (pictot_right != 1)
            {
              res2_right <- fpeaksFlat(spec_R,
                                       amp = amp,
                                       freq = freq,
                                       plot = plotpic,
                                       f)
              npic_right <- nrow(res2_right)
            } else{
              npic_right <- 1
            }
          } else{
            npic_right <- 0
          }
          Table_right <- cbind(Table_right, npic_right)
        }
      }
      
      
      Table_left <- as.data.frame(Table_left)
      Table_right <- as.data.frame(Table_right)
      if (stereo == FALSE)
      {
        Table_right <- Table_left
        Table_right[] <- NA
      }
      # modification tdelattre pour ne retourner que le channel left des wildlife acoustics mini bat
      # /!\ attention, c'est le channel right pour les mini "classic"
      #Result<-list(Table_left,Table_right)
      #names(Result)<-c("Mono_left","Mono_right")
      Result = Table_left
      return(Result)
    } else{
      #renvoie un dataframe remplit de NA pour ne pas faire planter la boucle.
      print("file too short, skipping and returning NA")
      #return(NA)
      return(data.frame(list(NA, NA, NA, NA, NA)))
    }
  }