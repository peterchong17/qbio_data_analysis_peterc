# exercise 1.1

#1.11
is.na(attenu)
#rows 79, 81, 94, 96, 99, 107, 108,114, 116, 118, 123, 126, 128

#1.12 
attenu_cleaned=na.omit(attenu)

#1.13
head(attenu_cleaned)
dim(attenu_cleaned)


# exercise 1.2

#1.21
Theoph_2 = Theoph
Theoph_2

#1.22
median(Theoph_2$Dose)

#1.23 
Theoph_2$Dose_Class=ifelse(Theoph_2$Dose>=median(Theoph_2$Dose), "high", "low")
Theoph_2

#1.24
print(head(Theoph_2))
print(dim(Theoph_2))


# exercise 1.3

#1.31
starbucks=read.csv("/Users/peterchong/Documents/qbio_data_analysis_peterc/week3_R/week3_homework/starbucks.csv")
starbucks

#1.32 
starbucks_bool=is.na(starbucks)
x=rowSums(starbucks_bool)
is_row_empty=x
is_row_empty[x==0] = TRUE
is_row_empty[x!=0]= FALSE
is_row_empty

if(length(is_row_empty)==nrow(starbucks)){
  print("yes")
}

starbucks_bool
starbucks_cleaned <- na.omit(starbucks)
starbucks_cleaned

#1.33
plot(x=starbucks_cleaned$Carb, y=starbucks_cleaned$Calories,
     main="Carbs vs Calories",
     xlab="Carbohydrates (g)",
     ylab="Calories")
#note that there is a positive correlation between carbohydrates and calories

#1.34
max_cal=max(starbucks_cleaned$Calories)
max_cal
l=which(starbucks_cleaned$Calories==max_cal)
l
is_highest_cal=starbucks_cleaned[y,]
is_highest_cal$Drink

#1.35
starbucks_cleaned$is_highest_fat=max(starbucks_cleaned$Fat)==starbucks_cleaned$Fat
plot( x = starbucks_cleaned$Carb, y = starbucks_cleaned$Calories,
     main="Carbs vs Calories",
     xlab="Carbohydrates (g)",
     ylab="Calories", col=factor(starbucks_cleaned$is_highest_fat))


# exercise 1.4

#1.41
baseball=read.csv("/Users/peterchong/Documents/qbio_data_analysis_peterc/week3_R/week3_homework/Batting.csv")

#1.42
z=which(baseball$HR>=3)
paste("There are ", length(z), " number of players that scored 3+ homeruns in a given year.")

#1.43
bball_graph=function(df){
  plot(x=df$yearID, y=df$HR,
     xlab="Year", ylab="# of Homeruns")
}
bball_graph(baseball)

#1.44
laangels=baseball[baseball$teamID=="LAA",]
bball_graph(laangels)

#1.45
atl_pit=baseball[baseball$teamID=="ATL"|baseball$teamID=="PIT",]
atl_pit
plot(x=atl_pit$yearID, y=atl_pit$HR,
     xlab="Year", ylab="# of Homeruns", col=factor(atl_pit$teamID))


# exercise 1.5

#1.51
easy_plot=function(a,b,color_data){
  m=median(color_data)
  levels=ifelse(color_data<m, "low", "high")
  levels=factor(levels)
  print(m)
#1.52
  plot(x=a, y=b,col=levels,
       pch=20)
#1.53
  cor.test(a,b)
}

#1.54
x=c(1,2,3,4,5,0)
y=c(5,4,3,2,1,8)
easy_plot(x, y, x)

easy_plot(starbucks_cleaned$Fat, starbucks_cleaned$Fiber, starbucks_cleaned$Sodium)
easy_plot(baseball$HR, baseball$G, baseball$SO)


# exercise 2.1
#The data set describes various quantitative measurements of several species of flowers' parts.
library(datasets)
data(iris)
summary(iris)
#It contains 150 observations of 5 variables.

# exercise 2.2
#Sepal length, sepal width, petal length, and petal width are all continuous variables. The species variable is categorical. The continuous variables are numeric, while the categorical is a string.

# exercise 2.3
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Petal.Width)
#I noticed that the two sepal variables are similar in shape, while the petal variables are also similar in shape. The petal histograms have an odd hole near the left of the graph.

# exercise 2.4
med_sep_width=median(iris$Sepal.Width)
iris_copy=iris
sepal_comparison=ifelse(iris_copy$Sepal.Width>med_sep_width, "wide", "narrow")
iris_copy$sep_comp=sepal_comparison
boxplot(iris_copy$Sepal.Width~iris_copy$sep_comp, ylab="Sepal Width")

# execise2.5
#The black species looks most unique, while the red and green species seem to overlap some and are therefore somewhat similar.
?pairs
args(pairs)
pairs(iris_copy[,1:4], col=iris_copy[,5])


