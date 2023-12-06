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
BCI2 <- BCI %>%
  select(date, Tair, wind.speed, SWdn, LWdn, RH, GPP) %>%
  rename(datetime = date) %>%
  mutate(RH = RH / 100,
         Tair = Tair + 273,
         date_str = lubridate::format_ISO8601(datetime, use_z_format=FALSE)) %>%
  select(date_str, everything(), -datetime) %>%
  drop_na()

MR2 <- MR %>%
  select(date, NEE, Tair, RH, Pressure, WS, SW_IN, SW_OUT, LW_IN, LW_OUT) %>%
  rename(datetime=date) %>%
  mutate(RH = RH / 100,
         Tair = Tair + 273,
         date_str = lubridate::format_ISO8601(datetime, use_z_format=FALSE)) %>%
  select(date_str, everything(), -datetime) %>%
  drop_na()

# Standardize names
names(BCI2) <- c(
  "date (:DateTime)", "Tair (K)", "wind_speed (m / s)", "SWdn (W / m^2)",
  "LWdn (W / m^2)", "RH", "GPP (μmol / m^2 / s)"
)

names(MR2) <- c(
  "date (:DateTime)", "GPP (μmol / m^2 / s)", "Tair (K)", "RH", "Pressure (kPa)",
  "wind_speed (m / s)", "SWdn (W / m^2)", "SWup (W / m^2)", "LWdn (W / m^2)",
  "LWup (W / m^2)"
)

# Write output
write_csv(BCI2, "data/BCI.csv")
write_csv(MR2, "data/MR.csv")

# Delete rdata files
system("bash -c 'rm -f data/BCI.Rdata'")
system("bash -c 'rm -f data/MR.Rdata'")
