library(tidyverse)
library(ggpubr)

as_tibble(iris)  # This is your data set for these problems
View(iris)  # View it in RStudio

### dplyr problems ###

# A. Complete the filter so that it only returns the "virginica" species
iris %>%
  filter(Species == "virginica")

# B. Use group_by and summarise to get the mean sepal length for each species
iris %>%
  group_by(Species) %>%
  summarise(sepal_mean = mean(Sepal.Length))

# C. Get the ratio of sepal length to sepal width for all samples
iris %>%
  mutate(sepal_ratio = Sepal.Length/Sepal.Width)

# D. Get the average ratio of petal length to petal width for each species
iris %>%
  mutate(petal_ratio = Petal.Length/Petal.Width) %>%
  group_by(Species) %>%
  summarise(petal_ratio_avg = mean(petal_ratio))

### ggplot problems ###

# A. Make a scatter plot comparing sepal length to petal length
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = Petal.Length)) +
  geom_point()

# B. Color the plot by the Species 
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point()

# C. Change to theme_classic(), set the base font size to 16, 
# set the x label to "Sepal Length (cm)", and set the y label to "Petal Length (cm)",
# set the title to "Sepal vs Petal Length in Iris Flowers"
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() +
  ylab("Petal Length (cm)") +
  xlab("Sepal Length (cm)") +
  labs(title = "Sepal vs Petal Length in Iris Flowers") +
  theme_classic(base_size = 16)

# D. Plot the linear correlation between Sepal and Petal length stat_smooth 
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() +
  ylab("Petal Length (cm)") +
  xlab("Sepal Length (cm)") +
  labs(title = "Sepal vs Petal Length in Iris Flowers") +
  theme_classic(base_size = 16) +
  stat_smooth(method = "lm")

# E. Add the pearson correlation coefficient using ggpubr's stat_cor()
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() +
  ylab("Petal Length (cm)") +
  xlab("Sepal Length (cm)") +
  labs(title = "Sepal vs Petal Length in Iris Flowers") +
  theme_classic(base_size = 16) +
  stat_smooth(method = "lm") +
  stat_cor()

# F. Save the plot to a pdf file
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() +
  ylab("Petal Length (cm)") +
  xlab("Sepal Length (cm)") +
  labs(title = "Sepal vs Petal Length in Iris Flowers") +
  theme_classic(base_size = 16) +
  stat_smooth(method = "lm") +
  stat_cor() +
  ggsave(filename = "my_figure.pdf")


### Challenge Problems ###

# Answer the following using dplyr and ggplot2 (with significance testing via ggpubr):

# 1. Do 6 cylinder cars have better MPG than 8 cylinder cars?
mtcars %>%
  filter(cyl == 6 | cyl == 8) %>%
  ggplot(mapping = aes(x = factor(cyl), y = mpg, fill = factor(cyl))) +
  geom_boxplot() +
  stat_compare_means(method = "t.test", label = "p.signif", 
                     comparisons = list(c("6","8"))) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none") +
  ylab("Miles per Gallon (MPG)") +
  xlab("Number of Cylinders") +
  labs(title = "Number of Cylinders and MPG")


# 2. Is there a correlation between urban population and Murder rate?
USArrests %>%
  ggplot(mapping = aes(x = UrbanPop, y = Murder)) + 
  geom_point() + 
  ylab("Murder rate (per 100k)") +
  xlab("Percentage of urban population") +
  labs(title = "Urban pop and per capita murder rate") +
  theme_classic() + 
  stat_smooth(method = "lm") +
  stat_cor()

# 3. Did the treatments make the plants grow more or less?
PlantGrowth %>%
  group_by(group) %>%
  ggplot(mapping = aes(x = factor(group, labels = c("Control", "Treatment 1", "Treatment 2")),
                       y = weight, fill = group)) +
  geom_boxplot() + 
  xlab(NULL) +
  ylab("Plant weight") +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  labs(title = "Plant Growth With Treatment") + 
  stat_compare_means(method = "t.test", label = "p.signif", 
                     comparisons = list(c("Control", "Treatment 1"), c("Control", "Treatment 2"))) 


# 4. Does Vitamin C supplementation improve tooth growth? Is OJ better than Vitamin C at dose of 2?
ToothGrowth %>%
  ggplot(aes(x = factor(supp), y = len, fill = supp)) +
  geom_boxplot() +
  facet_wrap(~ paste0(dose, " mg/day")) +
  stat_compare_means(method = "t.test",
                     label = "p.signif", 
                     comparisons = list(c("OJ", "VC"))) +
  xlab(NULL) +
  ylab("Tooth Length") +
  labs(title = "Tooth Growth With Treatment") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none") 




