# from https://www.youtube.com/watch?v=66z_MRwtFJM
source("DNAmAge_calc_full.R")
plot(all_data$age, all_data$Horvath, main="Scatterplot")
cor(all_data$age, all_data$Horvath) # 0.2939021

men1 = 'MEN1'
age_men1 = all_data$age[all_data$Sample_Group==men1]
age_men1_horvath = all_data$Horvath[all_data$Sample_Group==men1]
plot(age_men1, age_men1_horvath, main="Scatterplot")

cor(age_men1, age_men1_horvath) # 0.58267
mod_men1_horvath <- lm(age_men1_horvath ~ age_men1)
abline(mod_men1_horvath, col=2, lwd=3)
summary(mod_men1_horvath)
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -27.644  -9.672  -2.686   9.323  52.793 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  18.2026     9.7212   1.872   0.0685 .  
# age_men1      0.8155     0.1798   4.534 5.15e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 14.89 on 40 degrees of freedom
# Multiple R-squared:  0.3395,	Adjusted R-squared:  0.323 
# F-statistic: 20.56 on 1 and 40 DF,  p-value: 5.151e-05
             
sporadic = 'Sporadic'
age_sporadic = all_data$age[all_data$Sample_Group==sporadic]
age_sporadic_horvath = all_data$Horvath[all_data$Sample_Group==sporadic]
plot(age_sporadic, age_sporadic_horvath, main="Scatterplot")
cor(age_sporadic, age_sporadic_horvath) # -0.00519049

mod_sporadic_horvath <- lm(age_sporadic_horvath ~ age_sporadic)
abline(mod_sporadic_horvath, col=2, lwd=3)
summary(mod_sporadic_horvath)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -29.795 -12.446  -3.830   9.127  41.510 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  63.999379  10.422127   6.141 2.49e-07 ***
#   age_sporadic -0.006178   0.183665  -0.034    0.973    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 18.37 on 42 degrees of freedom
# Multiple R-squared:  2.694e-05,	Adjusted R-squared:  -0.02378 
# F-statistic: 0.001132 on 1 and 42 DF,  p-value: 0.9733

vhl = 'VHL'
age_vhl = all_data$age[all_data$Sample_Group==vhl]
age_vhl_horvath = all_data$Horvath[all_data$Sample_Group==vhl]
plot(age_vhl, age_vhl_horvath, main="Scatterplot")
cor(age_vhl, age_vhl_horvath) # 0.8753893

mod_vhl_horvath <- lm(age_vhl_horvath ~ age_vhl)
abline(mod_vhl_horvath, col=2, lwd=3)
summary(mod_vhl_horvath)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -8.826 -4.233 -2.643  6.386 11.476 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -16.5057    14.9785  -1.102 0.302526    
# age_vhl       1.6451     0.3212   5.122 0.000905 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 7.775 on 8 degrees of freedom
# Multiple R-squared:  0.7663,	Adjusted R-squared:  0.7371 
# F-statistic: 26.23 on 1 and 8 DF,  p-value: 0.0009052
