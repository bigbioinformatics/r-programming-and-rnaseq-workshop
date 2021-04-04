#using tidyverse to analyze qPCR Data
#Victoria Alers
#Bioinformatics Bootcamp


####Step 1: open all packages needed at begining######################

library(tidyverse)


#####Step 2: read in files############################################

#for tab delineated use:
res <- read_delim( "06202020.txt", "\t", 
                   escape_double = FALSE,
                   trim_ws = TRUE, skip = 9  )


# #for excel use:
# sams_data <- read_xlsx("06202020.xlsx",
#   #if the data is on a sheet other than the first of the workbook
#   sheet = 2,
#   #replaces all empty cells with the text: undetermined
#   na = "Undetermined",
#   #skip rows from 0-7; starts table at row 8
#   skip = 7,
#   #use only if you don`t want to import the whole sheet
#   n_max = 37)


#############Step 3: Tidy data#########################################

#without this all the arguments after will not save to the environment
tidyres <- res %>%
  #selects only the data you want
  select(Sample = "Sample Name", Primer = "Detector Name", Ct = "Ct") %>%
  #get rid of empty cells/rows
  drop_na() %>%
  #keeps rows that do not have Detector Name in Primer column
  #and Ct in the Ct column
  filter(Primer!="Detector Name", Ct!="Ct") %>%
  filter(Ct!="Undetermined") %>%
  #when read in because column Ct had characters 'undetermined'
  #we need to change the data type to use arithmetic
  #instead of making a new column, we use mutate_at
  mutate_at("Ct",as.numeric) 
#glimpse allows us to see how the data is saved in the environment
glimpse(tidyres)



###################################################################
#################### Hands on activity one ########################
###################################################################

## Import the data & preliminary cleaning ##

# Preliminary data -- primer key. Maps plate rows to primer info
primer_key <- data.frame(row = c("A", "B", "C", "D", "E", "F", "G", "H"),
                         primers = c(rep("HK", 4), rep("test", 4)))
treatment_key <- data.frame(column = as.character(1:12),
                            treatment = c(rep("RNAi1", 3),
                                          rep("RNAi2", 3),
                                          rep("Control", 3),
                                          rep(NA, 3)))

# Import the qPCR data
qPCR_data <- read_csv("https://raw.githubusercontent.com/liz-is/qpcr-analysis-with-r/gh-pages/_episodes_rmd/data/qpcr_data.csv",
                      skip = 27)

# Merge the primer_key and qPCR_data based on "row" column and...
# Also merge the treatment_key and qPCR_data based on "column" column
qPCR_data <- qPCR_data %>%
  mutate(row = substr(Well, 1, 1)) %>%
  mutate(column = substr(Well, 2, 2)) %>%
  left_join(y = primer_key, by = "row") %>%
  left_join(y = treatment_key, by = "column")

# Take a peak at the plate using ggplot2
ggplot(qPCR_data, aes(x = factor(column), y = row, 
                      fill = primers, label = `Sample Name`)) +
  geom_tile(colour = "black") +
  geom_text()

# Add a title: "My qPCR plate"
ggplot(qPCR_data, aes(x = factor(column), y = row, fill = primers, label = `Sample Name`)) +
  geom_tile(colour = "black") +
  geom_text() + 
  labs(title = "My qPCR plate")

# Change the x axis title to "Column" and y axis title to "Row"
ggplot(qPCR_data, aes(x = factor(column), y = row, fill = primers, label = `Sample Name`)) +
  geom_tile(colour = "black") +
  geom_text() + 
  labs(title = "My qPCR plate") + 
  ylab("Row") + 
  xlab("Column")

# Change the theme to theme_bw
ggplot(qPCR_data, aes(x = factor(column), y = row, fill = primers, label = `Sample Name`)) +
  geom_tile(colour = "black") +
  geom_text() + 
  labs(title = "My qPCR plate") + 
  ylab("Row") + 
  xlab("Column") +
  theme_bw()

# Color by "treatment" instead of "primers"
ggplot(qPCR_data, aes(x = factor(column), y = row, fill = treatment, label = `Sample Name`)) +
  geom_tile(colour = "black") +
  geom_text() + 
  labs(title = "My qPCR plate") + 
  ylab("Row") + 
  xlab("Column") +
  theme_bw()

## Tidy the data: ##

# 1. Select "Sample Name", "primers", "treatment", and "Ct" columns
# -- The new names of these columns should be "Sample", "Primer", "Treatment", and "Ct"
qPCR_data %>%
  _____(Sample = `Sample Name`, Primer = ______,
         ______ = ______, _____ = _____)

# 2. Remove all the rows with NA values in them
qPCR_data %>%
  _____(Sample = `Sample Name`, Primer = ______,
        ______ = ______, _____ = _____) %>%
  ______()

# 3. Remove all the rows with "Undetermined" as the Ct value
qPCR_data %>%
  _____(Sample = `Sample Name`, Primer = ______,
        ______ = ______, _____ = _____) %>%
  ______() %>%
  _______(Ct != "Undetermined")

# 4. Remove all the rows with "NTC" as the Sample
qPCR_data %>%
  _____(Sample = `Sample Name`, Primer = ______,
        ______ = ______, _____ = _____) %>%
  ______() %>%
  _______(Ct != "Undetermined") %>%
  ____________

# 5. Force Ct column to be a numeric value
qPCR_data %>%
  _____(Sample = `Sample Name`, Primer = ______,
        ______ = ______, _____ = _____) %>%
  ______() %>%
  _______(Ct != "Undetermined") %>%
  ____________ %>%
  _________("Ct",as.numeric) 

# 6. Assign your tidy data to a new tibble, "tidyPCR"
tidyPCR <- qPCR_data %>%
  _____(Sample = `Sample Name`, Primer = ______,
        ______ = ______, _____ = _____) %>%
  ______() %>%
  _______(Ct != "Undetermined") %>%
  ____________ %>%
  _________("Ct",as.numeric) 
tidyPCR

# 7. Glimpse tidyPCR
glimpse(tidyPCR)

###################################################################
###################################################################


#############Step 4: Transform data #######################################


# lets get the averages of each sample's house keeping gene (primer control)
HKgene <- tidyres %>%
  filter(Primer == "GAPDH") %>%
  group_by(Sample) %>%
  summarize(HKg = mean(Ct))


# Calculate ΔCt and expression 
dCt <- tidyres %>%
  filter(Primer != "GAPDH")%>%
  left_join(HKgene, by = "Sample") %>%
  mutate(dCt = Ct - HKg)


# get expression averages and descriptive statistics
# get the average of dCt
avgdCt <- dCt %>% 
  group_by(Sample, Primer) %>%
  summarise( avg_dCt = mean(dCt))

calcRes <- dCt %>%
  right_join(avgdCt, by= c("Sample", "Primer")) %>%
  group_by(Sample, Primer) %>%
  mutate(sdev = sd(dCt))

#############Step 5: Visualize your data ####################################


dCtplot <- ggplot(calcRes, aes(x = Sample, y = avg_dCt, fill = Primer)) + 
  geom_col( position = "dodge") +
  geom_errorbar(aes(ymin=avg_dCt-sdev, 
                    ymax=avg_dCt+sdev), 
                width=.2,
                position=position_dodge(.9))
dCtplot + scale_fill_viridis_d(option = "inferno")


############# Step 4: Transform data ddCt #########################################

# lets get the averages of each sample's house keeping gene dCt (primer control)
CTRdCt <- calcRes %>%
  filter(Sample == "U2OS_siscr") %>%
  group_by(Primer) %>%
  summarize(CTRdCt = mean(dCt))


# Calculate ddCt and expression 
ddCt <- calcRes %>%
  filter(Sample != "U2OS_siscr")%>%
  left_join(CTRdCt, by = "Primer") %>%
  mutate(ddCt = dCt - CTRdCt) 


# get ddCt averages and descriptive statistics
avgddCt <- ddCt %>% 
  group_by(Sample, Primer) %>%
  summarise( avg_ddCt = mean(ddCt))

calcRes_ddCt <- ddCt %>%
  right_join(avgddCt, by= c("Sample", "Primer")) %>%
  group_by(Sample, Primer) %>%
  mutate(sdev = sd(ddCt))

#############Step 5: Visualize your data ####################################

ddCtplot <- ggplot(calcRes_ddCt, aes(x = Sample, y = avg_ddCt, fill = Primer)) + 
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin=avg_ddCt-sdev, 
                    ymax=avg_ddCt+sdev), 
                width=.2,
                position=position_dodge(.9))
ddCtplot + scale_fill_viridis_d(option = "inferno")



###################################################################
#################### Hands on activity two ########################
###################################################################

# Recap: Assign your tidy data to a new tibble, "tidyPCR"
tidyPCR <- qPCR_data %>%
  _____(Sample = `Sample Name`, Primer = ______,
        ______ = ______, _____ = _____) %>%
  ______() %>%
  _______(Ct != "Undetermined") %>%
  ____________ %>%
  _________("Ct",as.numeric) 
tidyPCR

## Calculate dCt on tidyPCR ##

# Get the average CT for house keeping gene
HKgene <- tidyPCR %>%
  ________(Primer == "HK") %>%
  group_by(__________) %>%
  ________(______ = mean(Ct))

# Calculate ΔCt 
dCt <- tidyPCR %>%
  ______(Primer != "HK")%>%
  ______(y = HKgene, by = "Sample") %>%
  ______(dCt = Ct - HKg)

# get expression averages and descriptive statistics
# get the average of dCt
avgdCt <- dCt %>% 
  ____________ %>%
  __________(avg_dCt = mean(dCt))

calcRes <- dCt %>%
  right_join(avgdCt, by= c("Sample", "Primer", "treatment")) %>%
  __________________ %>%
  ___________(sdev = sd(dCt))

## Plot the dCt results ##

dCtplot <- _________(calcRes, aes(x = _____, y = ______, fill = Treatment)) + 
  _______( position = "dodge") +
  geom_errorbar(aes(ymin=avg_dCt-sdev, 
                    ymax=avg_dCt+sdev), 
                width=.2,
                position=position_dodge(.9)) +
  scale_fill_viridis_d(option = "inferno")
dCtplot

## Calculate ddCt and plot it ##

# lets get the averages of each sample's house keeping gene dCt (primer control)
CTRdCt <- calcRes %>%
  ________(Treatment == "Control") %>%
  pull(dCt) %>%
  mean()


# Calculate ddCt 
ddCt <- calcRes %>%
  ________(Treatment != "Control")%>%
  _________(ddCt = dCt - CTRdCt) 


# get ddCt averages and descriptive statistics
avgddCt <- ddCt %>% 
  _________(Sample) %>%
  __________(avg_ddCt = ________(_______))

calcRes_ddCt <- ddCt %>%
  __________(avgddCt, by= c("Sample")) %>%
  ________(________) %>%
  ________(sdev = sd(________))


## Plot the ddCt ##




###################################################################
###################################################################










