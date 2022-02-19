# from https://www.youtube.com/watch?v=66z_MRwtFJM
source("DNAmAge_calc_full.R")
png(file="/Users/michaelk/Documents/meth/pics/linear_regressions/linear_regression_all_horvath_age.png",
    width=800, height=600)
plot(all_data$age, all_data$Horvath, main="Multiple R-squared:  0.08638")
cor(all_data$age, all_data$Horvath) # 0.2939021
mod_all_horvath <- lm(all_data$Horvath ~ all_data$age)
abline(mod_all_horvath, col=2, lwd=3)
dev.off()
summary(mod_all_horvath)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -37.83 -10.18  -2.92   6.89  57.80 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   42.0630     6.9247   6.074 2.63e-08 ***
#   all_data$age   0.3781     0.1268   2.981  0.00366 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 17.05 on 94 degrees of freedom
# Multiple R-squared:  0.08638,	Adjusted R-squared:  0.07666 
# F-statistic: 8.887 on 1 and 94 DF,  p-value: 0.003656

men1 = 'MEN1'
age_men1 = all_data$age[all_data$Sample_Group==men1]
age_men1_horvath = all_data$Horvath[all_data$Sample_Group==men1]
png(file="/Users/michaelk/Documents/meth/pics/linear_regressions/linear_regression_men1_horvath_age.png",
    width=800, height=600)
plot(age_men1, age_men1_horvath, main="Multiple R-squared:  0.3395")

cor(age_men1, age_men1_horvath) # 0.58267
mod_men1_horvath <- lm(age_men1_horvath ~ age_men1)
abline(mod_men1_horvath, col=2, lwd=3)
dev.off()
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
png(file="/Users/michaelk/Documents/meth/pics/linear_regressions/linear_regression_sporadic_horvath_age.png",
    width=800, height=600)
plot(age_sporadic, age_sporadic_horvath, main="Multiple R-squared:  2.694e-05")
cor(age_sporadic, age_sporadic_horvath) # -0.00519049

mod_sporadic_horvath <- lm(age_sporadic_horvath ~ age_sporadic)
abline(mod_sporadic_horvath, col=2, lwd=3)
dev.off()
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

png(file="/Users/michaelk/Documents/meth/pics/linear_regressions/linear_regression_vhl_horvath_age.png",
    width=800, height=600)
plot(age_vhl, age_vhl_horvath, main="Multiple R-squared:  0.7663")
cor(age_vhl, age_vhl_horvath) # 0.8753893

mod_vhl_horvath <- lm(age_vhl_horvath ~ age_vhl)
abline(mod_vhl_horvath, col=2, lwd=3)
dev.off()
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

png(file="/Users/michaelk/Documents/meth/pics/linear_regressions/linear_regression_all_levine_age.png",
    width=800, height=600)
plot(all_data$age, all_data$Levine, main="Multiple R-squared:  0.14278")
cor(all_data$age, all_data$Levine) # 0.3777282
mod_all_levine <- lm(all_data$Levine ~ all_data$age)
abline(mod_all_levine, col=2, lwd=3)
dev.off()
summary(mod_all_levine)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -45.465 -10.499   0.564  11.924  37.445 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    7.0417     7.0640   0.997 0.321400    
# all_data$age   0.5118     0.1294   3.955 0.000148 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 17.39 on 94 degrees of freedom
# Multiple R-squared:  0.1427,	Adjusted R-squared:  0.1336 
# F-statistic: 15.64 on 1 and 94 DF,  p-value: 0.0001482

age_men1_levine = all_data$Levine[all_data$Sample_Group==men1]
png(file="/Users/michaelk/Documents/meth/pics/linear_regressions/linear_regression_men1_levine_age.png",
width=800, height=600)
plot(age_men1, age_men1_levine, main="Multiple R-squared:  0.2549")

cor(age_men1, age_men1_levine) # 0.5048464
mod_men1_levine <- lm(age_men1_levine ~ age_men1)
abline(mod_men1_levine, col=2, lwd=3)
dev.off()
summary(mod_men1_levine)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -37.483  -6.075   2.428  11.660  27.810 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -7.5513    10.1621  -0.743 0.461778    
# age_men1      0.6954     0.1880   3.699 0.000651 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 15.56 on 40 degrees of freedom
# Multiple R-squared:  0.2549,	Adjusted R-squared:  0.2362 
# F-statistic: 13.68 on 1 and 40 DF,  p-value: 0.0006506

age_sporadic_levine = all_data$Levine[all_data$Sample_Group==sporadic]
png(file="/Users/michaelk/Documents/meth/pics/linear_regressions/linear_regression_sporadic_levine_age.png",
    width=800, height=600)
plot(age_sporadic, age_sporadic_levine, main="Multiple R-squared:  0.06854")

cor(age_sporadic, age_sporadic_levine) # 0.2617971
mod_sporadic_levine <- lm(age_sporadic_levine ~ age_sporadic)
abline(mod_sporadic_levine, col=2, lwd=3)
dev.off()
summary(mod_sporadic_levine)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -49.066 -12.805   2.621  13.491  33.425 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   20.4206    10.9246   1.869   0.0686 .
# age_sporadic   0.3384     0.1925   1.758   0.0860 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 19.26 on 42 degrees of freedom
# Multiple R-squared:  0.06854,	Adjusted R-squared:  0.04636 
# F-statistic:  3.09 on 1 and 42 DF,  p-value: 0.08604
# 

age_vhl_levine = all_data$Levine[all_data$Sample_Group==vhl]
png(file="/Users/michaelk/Documents/meth/pics/linear_regressions/linear_regression_vhl_levine_age.png",
    width=800, height=600)
plot(age_vhl, age_vhl_levine, main="Multiple R-squared:  0.5309")

cor(age_vhl, age_vhl_levine) # 0.7286316
mod_vhl_levine <- lm(age_vhl_levine ~ age_vhl)
abline(mod_vhl_levine, col=2, lwd=3)
dev.off()
summary(mod_vhl_levine)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8.9723 -7.7871 -0.0491  4.6256 12.7620 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -15.9271    16.8924  -0.943   0.3734  
# age_vhl       1.0900     0.3622   3.009   0.0168 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.769 on 8 degrees of freedom
# Multiple R-squared:  0.5309,	Adjusted R-squared:  0.4723 
# F-statistic: 9.054 on 1 and 8 DF,  p-value: 0.01684
