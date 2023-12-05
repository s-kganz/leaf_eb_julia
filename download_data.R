library(tidyverse)
library(lubridate)

BCI_link <- "https://zenodo.org/records/6862565/files/BCI.RData"
MR_link  <- "https://zenodo.org/records/6862565/files/MR.RData"

download.file(BCI_link, "data/BCI.Rdata")
download.file(MR_link, "data/MR.Rdata")

load("data/BCI.Rdata")
load("data/MR.Rdata")

# Select important columns. These aren't really documented on Zenodo so I'm
# guessing on what some of these mean.
BCI <- BCI %>%
  select(date, DOY, Tcan.full, Tcan.roi.evergreen, Tcan.roi.deciduous,
         Tair, wind.speed, SWdn, LWdn, RH, GPP)

MR <- MR %>%
  select(DateTime, NEE, Tair, RH, Pressure, WS, SW_IN, SW_OUT, LW_IN, LW_OUT)

# Standardize names
names(BCI) <- c(
  "date (:Date)", "DOY (day)", "Tair (°C)", "wind.speed (m / s)", "SWdn (W / m^2 / s)",
  "LWdn (W m^-2 s^-1)", "RH (%)", "GPP (μmol / m^2 / s)"
)

names(MR) <- c(
  "date (:Date)", "GPP (μmol / m^2 / s)", "Tair (°C)", "RH (%)", "Pressure (kPa)",
  "wind.speed (m / s)", "SWdn (W / m^2 / s)", "SWup (W / m^2 / s)", "LWdn (W / m^2 / s)",
  "LWup (W / m^2 / s)"
)

# Write output
write_csv(BCI, "data/BCI.csv")
write_csv(MR, "data/MR.csv")

# Delete rdata files
system("bash -c 'rm -f data/BCI.Rdata'")
system("bash -c 'rm -f data/MR.Rdata'")
