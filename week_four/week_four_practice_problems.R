library(tidyverse)
library(ggpubr)

as_tibble(iris)  # This is your data set for these problems
View(iris)  # View it in RStudio

### dplyr problems ###

# NOTE: Don't forget, you there's a cheat sheet available on Box under "Resources"

# A. Complete the filter so that it only returns the "virginica" species
iris %>%
  filter(___)

# B. Use summarise to get the mean sepal length for each species
iris %>%
  group_by(____) %>%
  summarise(sepal_mean = mean(____))

# C. Get the ratio of sepal length to sepal width for all samples
iris %>%
  _____(sepal_ratio = _____/Sepal.Width)

# D. Get the average ratio of petal length to petal width for each species
iris %>%
  _____(petal_ratio = _____/Petal.Width) %>%
  _____(_____) %>%
  ______(petal_ratio_avg = _____(petal_ratio))

### ggplot problems ###

# A. Make a scatter plot comparing sepal length to petal length
iris %>%
  ggplot(mapping = ___(x = Sepal.Length, y = ____)) +
  geom______()

# B. Color the plot by the Species 
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = _____, _____ = _____)) +
  geom______()

# C. Change to theme_classic(), set the base font size to 16, 
# set the x label to "Sepal Length (cm)", and set the y label to "Petal Length (cm)",
# set the title to "Sepal vs Petal Length in Iris Flowers"
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = ______, _____ = _____)) +
  geom______() +
  _____("Petal Length (cm)") +
  _____(_______) +
  ______(title = "Sepal vs Petal Length in Iris Flowers") +
  theme_______(_______)

# D. Plot the linear correlation between Sepal and Petal length stat_smooth 
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = ______, _____ = _____)) +
  geom______() +
  _____("Petal Length (cm)") +
  _____(_______) +
  ______(title = "Sepal vs Petal Length in Iris Flowers") +
  theme_______(_______) +
  stat_smooth(method = _____)

# E. Add the pearson correlation coefficient using ggpubr's stat_cor()
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = ______, _____ = _____)) +
  geom______() +
  _____("Petal Length (cm)") +
  _____(_______) +
  ______(title = "Sepal vs Petal Length in Iris Flowers") +
  theme_______(_______) +
  stat_smooth(method = _____) +
  stat_cor()

# F. Save the plot to a pdf file
iris %>%
  ggplot(mapping = aes(x = Sepal.Length, y = ______, _____ = _____)) +
  geom______() +
  _____("Petal Length (cm)") +
  _____(_______) +
  ______(title = "Sepal vs Petal Length in Iris Flowers") +
  theme_______(_______) +
  stat_smooth(method = _____) +
  stat_cor() +
  ______(________ = "my_figure.pdf")


### Challenge Problems ###

# Answer the following using dplyr and ggplot2 (with significance testing via ggpubr):

# 1. Do 6 cylinder cars have better MPG than 8 cylinder cars?
View(mtcars)


# 2. Is there a correlation between urban population and Murder rate?
View(USArrests)


# 3. Did the treatment #1 or #2 make the plants grow more or less than control?
View(PlantGrowth)


# 4. Does Vitamin C supplementation improve tooth growth? Is OJ better than Vitamin C at dose of 2?
View(ToothGrowth)


