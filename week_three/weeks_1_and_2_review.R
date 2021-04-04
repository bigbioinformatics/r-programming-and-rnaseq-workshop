#### Weeks One and Two Review ####

### Introductory concepts ###

## Simple data structures

# Numerics
1
1E100
1.1


# Logicals
TRUE
FALSE


# Characters
"Hello world!"
"This is a string!"


## Simple data structure methods

# Class
class("Hello world!")
class(1)
class(TRUE)

# Structure
str("Hello world!")
str(1)
str(TRUE)


# Arithmetic
1 + 1  # add
1 - 1  # subtract
2 * 5  # multiply
10 / 5  # divide
2 ^ 3  # exponent
10 %% 3  # modulo
14 %% 3

# Logical operations
! TRUE  # "not" TRUE is FALSE
! FALSE # "not" FALSE is TRUE

went_to_store <- TRUE
bought_cheese <- TRUE

# Did you go to the store AND buy cheese?
(went_to_store & bought_cheese)

# Did you get up on time OR take the bus?
up_on_time <- FALSE
took_the_bus <- FALSE
(up_on_time | took_the_bus)

# Did you get (flowers AND potatoes OR bread)? 
flowers <- TRUE
potatoes <- FALSE
bread <- FALSE
((flowers & potatoes) | ! bread)

TRUE & FALSE  # "&" (and) is TRUE only when both sides are TRUE
TRUE | FALSE  # "|" (or) is TRUE when at least one side is TRUE

TRUE && FALSE && "Hello world!"
FALSE || FALSE || "Hello world!"



# Comparison (works with any data type)
(1 == 1)  # equivalence -- is "1" equal to "1"?
"Hello" == "World"
10 != 1  # non-equivalence
10 > 1  # greater-than
10 < 100 # less-than
2 <= 2  # less-than or equal to
5 >= 3  # greater-than or equal to


## Complex Data Structures

# Variables
a <- 1
b <- 2
a + b
# variables can have names
names(a) <- "Hello world!"
a


# Vectors
c(1, 2, 3)
c("Hello", "World", "!")
c(TRUE, FALSE, FALSE)
# You can make a numeric vector with shorthand :
num_vec <- 1:10
num_vec
# Vectors can have names
my_vec <- c(1, 2, 3)
names(my_vec) <- c("Hello", "World", "!")
my_vec
# Elements of a vector can be access by index
my_vec[1]  # Get index 1
my_vec[1:2]  # Get indices 1->2
# Elements of a vector can be accessed by name
my_vec["Hello"]  # Get the element with name "Hello"


# Data Frames
my_df <- data.frame(
  "Numerics" = c(1, 2, 3),
  "Characters" = c("Hello", "World", "!"),
  "Logicals" = c(TRUE, FALSE, FALSE)
)
my_df
# Data Frames can have row names and column names
rownames(my_df) <- c("Row_1", "Row_2", "Row_3")
colnames(my_df) <- c("Col_1", "Col_2", "Col_3")
my_df
# Data Frames can be viewed in RStudio
View(my_df)
# Data Frames can be accessed by row and column index
my_df[1,2]  # Row 1, column 2
my_df[2,1]  # Row 2, column 1
my_df[3,]  # Row 3, all columns
my_df[,2]  # all rows, column 2
# Data Frames can be accessed by row and column name
my_df["Row_1", 2]  # Row 1, column 2
my_df[1, "Col_3"]  # Row 1, column 3
my_df["Row_2", "Col_1"]  # Row 2, column 1
my_df["Row_2",]  # Row 2, all columns
# Data Frame columns can be accessed as a vector using $ notation
my_df$Col_1
my_df$Col_2[2]  # Column 2, row 2


# Lists
my_list <- list(1, "a", TRUE, c(1, 2, 3))
my_list
# Lists can have names
names(my_list) <- c("Elem_1", "Elem_2", "Elem_3", "Elem_4")
my_list
# Lists can be accessed by element index
my_list[[2]]  # This returns the data
# Lists can be accessed by name
my_list[["Elem_3"]]
# Lists can be accessed by name using $ notation
my_list$Elem_4


# Matrices
my_mat <- matrix(1:10, nrow = 2, ncol = 5, byrow = TRUE)
my_mat
# Matrices can have column and row names
colnames(my_mat) <- c("Col 1", "Col 2", "Col 3", "Col 4", "Col 5")
rownames(my_mat) <- c("Row 1", "Row 2")
my_mat
# Matrices can be accessed using indices or row/col names
my_mat[1, "Col 1"]  # Row 1, column 1
my_mat["Row 2", "Col 5"]  # Row 2, column 5
my_mat[2, 3]  # Row 2, column 3


# Factors 
my_fct <- factor(x = c(90, 80, 70, 60, 0))
my_fct
# Factors can be named
names(my_fct) <- c("A", "B", "C", "D", "F")
my_fct
# Factors can be accessed the same as vectors
my_fct[1]
my_fct["D"]
# Factors can be re-leveled
levels(my_fct)
levels(my_fct) <- c(90, 80, 70, 60, 0)
levels(my_fct)


### Intermediate concepts ###

## Complex Data Structure Methods

# Arithmetic
vec1 <- c(1, 10, 4, 5)
vec1 + 1  # Adds 1 to every element
vec2 <- c(10, 2, 3)
vec1 + vec2  # Pairwise addition
vec1 * vec2  # Pairwise multiplication


# Comparison
vec1 <- c(1, 10, 4)
vec1 > 1  # Element-wise greater-than 
vec2 <- c(10, 2, 3)
vec1 > vec2  # Pairwise less than
vec1 < vec2  # Pairwise greater than
vec1 == vec2  # Pairwise equivalence
vec1 != vec2  # Pairwise non equivalence
vec3 <- c(TRUE, FALSE, TRUE)
! vec3  # Negate every element
vec4 <- c(FALSE, TRUE, TRUE)
vec3 & vec4  # Requires both elements are TRUE
vec3 | vec4  # Only one element must be TRUE
vec5 <- c("Hello", "World", "!")
vec5 == "Hello"
vec6 <- c("Aloha", "World", "!")
vec5 == vec6

# The %in% operator
11 %in% 1:10
5 %in% 1:10
TRUE %in% c(TRUE, FALSE, FALSE)
"Hello" %in% c("Hello", "World", "!")
! "Aloha" %in% c("Hello", "World", "!") 

# The which() function
1:10 == 5
which(1:10 == 5)  # Return the index of all "TRUE"s in a logical vector

# The length() function
my_vec <- 1:5
length(my_vec)  # Gives the length of the vector

# The sample() function
my_vec <- 1:100
sample(my_vec, size = 10)  # Get 10 random values from my_vec

# max, min, mean, sum, median, etc function family
my_nums <- c(5, 100, 20, 4, 1.5, 3.333)
max(my_nums)  # Returns max value
min(my_nums)  # Returns min value
mean(my_nums)  # Returns mean value
median(my_nums)  # Returns median value
sum(my_nums)  # Returns sum of numeric vector
sd(my_nums)  # Returns standard deviation of numeric vector

## If logic

# If statements
i <- 4
if (i > 5) {
  print("i is greater than 5!")
}

# If...Else statements
i <- 4
if (i > 5) {
  print("i is greater than 5!")
} else {
  print("i is NOT greater than 5!")
}

# If...Else if...Else statements
i <- 0
if (i > 0) {
  print("i is positive!")
} else if (i == 0) {
  print("i is zero!")
} else {
  print("i is negative!")
}


## Loops

# For loops
days_of_week <- c("Monday", "Tuesday", "Wednesday", "Thursday", 
                  "Friday", "Saturday", "Sunday")
# Type I
for (day in days_of_week) {
  print(day)
}
# Type II
for (i in 1:length(days_of_week)) {
  day <- days_of_week[i]
  print(day)
}

# While loops
x <- 5
while (x > 0) {
  print(x)
  x <- x - 1
}

## Combining loops and if...Else
days_of_week <- c("Monday", "Tuesday", "Wednesday", "Thursday", 
                  "Friday", "Saturday", "Sunday")
for (day in days_of_week) {
  print(day)
  if (day == "Saturday") {
    print("Party time! :)")
  } else if (day == "Monday") {
    print("I HATE mondays! :( :(")
  } else {
    print("No parties today... :(")
  }
}

## Functions

# Functions take an input and perform an operation
my_funct <- function (number) {
  print("Your number is:")
  print(number)
}
my_funct(100)

# Functions may return an output
random_picker <- function(a_vector) {
  result <- sample(a_vector, size = 1)
  return(result)
}
random_picker(days_of_week)


## The apply family

# sapply
my_function <- function(day) {
  if (day == "Saturday") {
    return("Party time! :)")
  } else if (day == "Monday") {
    return("I HATE mondays! :( :(")
  } else {
    return("No parties today... :(")
  }
}
res_vec <- sapply(days_of_week, FUN = function(day) {
  if (day == "Saturday") {
    return("Party time! :)")
  } else if (day == "Monday") {
    return("I HATE mondays! :( :(")
  } else {
    return("No parties today... :(")
  }
})
res_vec

# lapply
res_lst <- lapply(days_of_week, FUN = function(day) {
  if (day == "Saturday") {
    return("Party time! :)")
  } else if (day == "Monday") {
    return("I HATE mondays! :( :(")
  } else {
    return("No parties today... :(")
  }
})
names(res_lst) <- days_of_week
res_lst

# apply (for dataframes/matrices)
my_mat <- matrix(1:10, nrow = 2, ncol = 5)
my_mat
# Apply across all matrix elements
res_mat <- apply(my_mat, MARGIN = 1:2, FUN = function(number) {
  if (number > 5) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})
res_mat
# Apply across columns
res_mat <- apply(my_mat, MARGIN = 2, function(column) {
  return(sum(column))
})
res_mat
# Apply across rows
res_mat <- apply(my_mat, MARGIN = 1, function(row) {
  return(sum(row))
})
res_mat













