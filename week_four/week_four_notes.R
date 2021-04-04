#### Week Three Lecture Notes #####

## Getting started with Tidyverse: Readr and Dplyr ##

# Download tidyverse and ggpubr
install.packages("tidyverse")
install.packages("ggpubr")

# Load the tidyverse package
library(tidyverse)

# Grab the titanic data
tidy_titanic <- read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv")
tidy_titanic

# Filter the titanic data for passengers > age 12
tidy_titanic %>%
  filter(Age > 12)

# Select the "Survived" and "Sex" columns
tidy_titanic %>%
  filter(Age > 12) %>%
  select(Survived, Sex)

# Group by Sex, then Survival
tidy_titanic %>%
  filter(Age > 12) %>%
  select(Survived, Sex) %>%
  group_by(Sex, Survived)

# Summarise
tidy_titanic %>%
  filter(Age > 12) %>%
  select(Survived, Sex) %>%
  group_by(Sex, Survived) %>%
  summarise(num_sex_survive = n())

# Or... just use tally
tidy_titanic %>%
  filter(Age > 12) %>%
  select(Survived, Sex) %>%
  group_by(Sex, Survived) %>%
  tally()

# Review... what does sum() do?
my_nums <- c(1, 5, 8)
sum(my_nums)  # sums
my_nums / sum(my_nums)  # proportions
(my_nums / sum(my_nums)) * 100  # percentages


# Mutate to change count to percentage of the sex surviving/dying
tidy_titanic %>%
  filter(Age > 12) %>%
  select(Survived, Sex) %>%
  group_by(Sex, Survived) %>%
  summarise(num_sex_survive = n()) %>%
  mutate(pct_sex_survive = num_sex_survive/sum(num_sex_survive) * 100)

# Keep only the percentage of those surviving
tidy_titanic %>%
  filter(Age > 12) %>%
  select(Survived, Sex) %>%
  group_by(Sex, Survived) %>%
  summarise(num_sex_survive = n()) %>%
  mutate(pct_sex_survive = num_sex_survive/sum(num_sex_survive) * 100) %>%
  filter(Survived == 1)

# Save the result as a variable
pct_surv <- tidy_titanic %>%
  filter(Age > 12) %>%
  select(Survived, Sex) %>%
  group_by(Sex, Survived) %>%
  summarise(num_sex_survive = n()) %>%
  mutate(pct_sex_survive = num_sex_survive/sum(num_sex_survive) * 100) %>%
  filter(Survived == 1)
pct_surv


## Hands-on activity One ##

View(mtcars)

# Filter mtcars for cars where "cyl" is 6 or more
mtcars %>%
  filter(cyl >= 6)

# Select the "disp" and "wt" and "mpg" columns from mtcars
mtcars %>%
  select(disp, wt, mpg)

# Group mtcars by number of cylinders, "cyl"
mtcars %>%
  group_by(cyl)

# Add up the number of cars in each cylinder group
mtcars %>%
  group_by(cyl) %>%
  summarize(num_cars = n())

# Calculate the median weight, "wt", in each cylinder group
mtcars %>%
  group_by(cyl) %>%
  summarize(median_wt = median(wt))

# Create a new column which contains the ratio of horse power, "hp", to "mpg"
# AND only keep the cars with "wt" greater than 3
mtcars %>%
  mutate(hp_ratio = hp / mpg) %>%
  filter(wt > 3)

# Group mtcars by "gear" and then find the average "hp" to "wt" ratio
mtcars %>%
  mutate(hp_ratio = hp / wt) %>%
  group_by(gear) %>%
  summarize(avg_hp_ratio = mean(hp_ratio))


## Getting started with Tidyverse: Readr and Dplyr ##

pct_surv <- tidy_titanic %>%
  filter(Age > 12) %>%
  select(Survived, Sex) %>%
  group_by(Sex, Survived) %>%
  summarise(num_sex_survive = n()) %>%
  mutate(pct_sex_survive = num_sex_survive/sum(num_sex_survive) * 100) %>%
  filter(Survived == 1)
pct_surv

# Set the data layer
ggplot(data = pct_surv)

# Set the aesthetic layer
ggplot(data = pct_surv, mapping = aes(x = Sex, y = num_sex_survive))

# Set the geometry layer
ggplot(data = pct_surv, mapping = aes(x = Sex, y = num_sex_survive)) + 
  geom_bar(stat = "identity")

# Change the y-axis aesthetic
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive)) + 
  geom_bar(stat = "identity")

# Change the y-axis aesthetic
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive)) + 
  geom_bar(stat = "identity")

# Add y-axis label
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving")

# Add title
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving") +
  labs(title = "Effect of Sex on Survival aboard the Titanic")

# Change the theme
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving") +
  labs(title = "Effect of Sex on Survival aboard the Titanic") +
  theme_classic()

# Color the bars by the Sex column
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive, fill = Sex)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving") +
  labs(title = "Effect of Sex on Survival aboard the Titanic") +
  theme_classic() 

# Remove the x-axis label
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive, fill = Sex)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving") +
  labs(title = "Effect of Sex on Survival aboard the Titanic") +
  theme_classic() +
  xlab(NULL)

# Remove the legend
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive, fill = Sex)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving") +
  labs(title = "Effect of Sex on Survival aboard the Titanic") +
  theme_classic() +
  xlab(NULL) +
  theme(legend.position = "none")

# Change the y-axis range and remove the expand
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive, fill = Sex)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving") +
  labs(title = "Effect of Sex on Survival aboard the Titanic") +
  theme_classic() +
  xlab(NULL) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

# Capitalize the x-axis labels
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive, fill = Sex)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving") +
  labs(title = "Effect of Sex on Survival aboard the Titanic") +
  theme_classic() +
  xlab(NULL) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Female", "Male"))

# Increase the font size
ggplot(data = pct_surv, mapping = aes(x = Sex, y = pct_sex_survive, fill = Sex)) + 
  geom_bar(stat = "identity") +
  ylab("Percentage Surviving") +
  labs(title = "Effect of Sex on Survival aboard the Titanic") +
  theme_classic(base_size = 15) +
  xlab(NULL) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Female", "Male"))

## Hands-on activity Two ##

View(iris)

# Group iris by "Species" 

iris %>%
  group_by(Species)

# Find the average "Sepal.length" to "Petal.length" ratio per species

iris %>%
  group_by(Species) %>%
  summarise(len_ratio = mean(Sepal.Length / Petal.Length))
  

# Pipe this into ggplot2. Set x axis to "Species", y to the ratio column

iris %>%
  group_by(Species) %>%
  summarise(len_ratio = mean(Sepal.Length / Petal.Length)) %>%
  ggplot(mapping = aes(x = Species, y = len_ratio))

# Add the barplot geometry

iris %>%
  group_by(Species) %>%
  summarise(len_ratio = mean(Sepal.Length / Petal.Length)) %>%
  ggplot(mapping = aes(x = Species, y = len_ratio)) +
  geom_bar(stat = "identity")

# Set y axis label to "Average of Sepal to Petal length (cm)"

iris %>%
  group_by(Species) %>%
  summarise(len_ratio = mean(Sepal.Length / Petal.Length)) %>%
  ggplot(mapping = aes(x = Species, y = len_ratio)) +
  geom_bar(stat = "identity") +
  ylab("Average of Sepal to Petal length (cm)")

# Color the plot by species


iris %>%
  group_by(Species) %>%
  summarise(len_ratio = mean(Sepal.Length / Petal.Length)) %>%
  ggplot(mapping = aes(x = Species, y = len_ratio, fill = Species)) +
  geom_bar(stat = "identity") +
  ylab("Average of Sepal to Petal length (cm)")


# Set the title to "Iris flower analysis"

iris %>%
  group_by(Species) %>%
  summarise(len_ratio = mean(Sepal.Length / Petal.Length)) %>%
  ggplot(mapping = aes(x = Species, y = len_ratio, fill = Species)) +
  geom_bar(stat = "identity") +
  ylab("Average of Sepal to Petal length (cm)") +
  labs(title = "Iris flower analysis")



## Getting started with Tidyverse: Readr and Dplyr Continued... ##

# Set x axis to Survived and y to Fare
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  ggplot(mapping = aes(x = Survived, y = Fare)) 

# Survived is being treated as a continuous variable...
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  select(Survived)

# Make it discrete with factor()
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  ggplot(mapping = aes(x = factor(Survived), y = Fare)) 

# Set factor labels to improve readability
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  ggplot(mapping = aes(x = factor(Survived, labels = c("Died", "Survived")), y = Fare)) 

# Add boxplot geometry layer
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  ggplot(mapping = aes(x = factor(Survived, labels = c("Died", "Survived")), y = Fare)) +
  geom_boxplot() 

# Log scale for y-axis
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  filter(Fare > 0) %>%
  ggplot(mapping = aes(x = factor(Survived, labels = c("Died", "Survived")), y = Fare)) +
  geom_boxplot() +
  scale_y_log10()

# Bracket passenger age into "Senior", "Adult" and "Child" groups
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  
  filter(is.numeric(Age) & Fare > 0) %>%
  mutate(Age_bracket = case_when(
    Age > 60 ~ "Senior",
    Age < 16 ~ "Child",
    TRUE ~ "Adult"  # This TRUE is equivalent to "else" in if...else
  ))

# Use facets to split the plot based on these brackets 
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  filter(is.numeric(Age) & Fare > 0) %>%
  mutate(Age_bracket = case_when(
    Age > 60 ~ "Senior",
    Age < 16 ~ "Child",
    TRUE ~ "Adult"  # This TRUE is equivalent to "else" in if...else
  )) %>%
  ggplot(mapping = aes(x = factor(Survived, labels = c("Died", "Survived")), y = Fare)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~ Age_bracket) 

# Load ggpubr
library(ggpubr)

# Add a t-test to our plot
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  filter(is.numeric(Age) & Fare > 0) %>%
  mutate(Age_bracket = case_when(
    Age > 60 ~ "Senior",
    Age < 16 ~ "Child",
    TRUE ~ "Adult"  # This TRUE is equivalent to "else" in if...else
  )) %>%
  ggplot(mapping = aes(x = factor(Survived, labels = c("Died", "Survived")), y = Fare)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~ Age_bracket) +
  stat_compare_means(method = "t.test")


# Make the label more reader-friendly 
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  filter(is.numeric(Age) & Fare > 0) %>%
  mutate(Age_bracket = case_when(
    Age > 60 ~ "Senior",
    Age < 16 ~ "Child",
    TRUE ~ "Adult"  # This TRUE is equivalent to "else" in if...else
  )) %>%
  ggplot(mapping = aes(x = factor(Survived, labels = c("Died", "Survived")), y = Fare)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~ Age_bracket) +
  stat_compare_means(method = "t.test", label = "p.signif")

# User "comparisons" to improve the readability even more
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  filter(is.numeric(Age) & Fare > 0) %>%
  mutate(Age_bracket = case_when(
    Age > 60 ~ "Senior",
    Age < 16 ~ "Child",
    TRUE ~ "Adult"  # This TRUE is equivalent to "else" in if...else
  )) %>%
  ggplot(mapping = aes(x = factor(Survived, labels = c("Died", "Survived")), y = Fare)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~ Age_bracket) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Died", "Survived")))


# Improve the overall appearance
read_csv("https://gist.githubusercontent.com/michhar/2dfd2de0d4f8727f873422c5d959fff5/raw/fa71405126017e6a37bea592440b4bee94bf7b9e/titanic.csv") %>%  # Read the data into a tibble
  filter(is.numeric(Age) & Fare > 0) %>%
  mutate(Age_bracket = case_when(
    Age > 60 ~ "Senior",
    Age < 16 ~ "Child",
    TRUE ~ "Adult"  # This TRUE is equivalent to "else" in if...else
  )) %>%
  ggplot(mapping = aes(x = factor(Survived, labels = c("Died", "Survived")), y = Fare, fill =  factor(Survived, labels = c("Died", "Survived")))) +
  geom_boxplot() +
  scale_y_log10(limits = c(2, 1000)) +
  facet_wrap(~ Age_bracket) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Died", "Survived"))) +
  xlab(NULL) +
  ylab("Ticket Price (log scale)") +
  theme_bw(base_size = 17) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("goldenrod", "skyblue")) +
  ggsave(filename = "my_figure.png", height = 5, width = 11)








