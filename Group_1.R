
# Group 1 -------------------------------------------------------------------
# Alec Santiago, Jay Mantuhac, Lynn Gao
# Stats 112


# Exploratory Data Analysis -----------------------------------------------

library(tidyverse)
library(nlme)
library(mgcv)
library(geepack)
library(glme)
library(lme4)
aids = read.csv("aids.csv")

glimpse(aids)

subjects = unique(aids$id)
subjects = max(subjects)

summary(aids)

aids %>% 
  summarize(mean_log_cd4 = mean(log_cd4, na.rm = TRUE), 
            stdev_log_cd4 = sd(log_cd4, na.rm = TRUE))
aids %>% 
  group_by(treatment) %>% 
  summarize(mean_log_cd4 = mean(log_cd4, na.rm = TRUE), 
            stdev_log_cd4 = sd(log_cd4, na.rm = TRUE))

# Distribution of Age
aids %>% 
  ggplot(aes(x = age)) + 
  geom_histogram(bins = 15, color = "black", fill = "gray")

# Distribution of people in each treatment group
table(aids$treatment)
round(prop.table(table(aids$treatment)), 2)

## Bivariate summaries 

# Boxplot based on treatment group
aids %>% 
  ggplot(aes(x = factor(treatment), y = log_cd4, 
             fill = factor(treatment))) +
  geom_boxplot()

# Boxplot based on gender
aids %>% 
  ggplot(aes(x = gender, y = log_cd4, fill = gender)) +
  geom_boxplot()

## Overall trends of response 

# Mean log_cd4 per week of each treatment 
aids %>% 
  ggplot(aes(x = week, y = log_cd4,  color =  factor(treatment))) +  
  labs(title = "Smoothed Trend of Each Treatment",
       y = "Log cd4", x = "Week") +
  geom_smooth(se = FALSE) +
  facet_wrap(~factor(treatment))

# Trends pf each subject per week separated by treatment
aids %>% 
  ggplot(aes(x = week, y = log_cd4, color = factor(treatment)))+
  geom_point(alpha = 0.3) +
  geom_smooth(se = FALSE) +
  coord_cartesian(ylim = c(2.2,3.3)) + 
  labs(title = "Smoothed Trend of Each Subject",
       y = " Log cd4", x = "Week")

aids %>% 
  group_by(treatment, week) %>% 
  summarize(mean_log_cd4 = mean(log_cd4, na.rm = TRUE)) %>% 
  ggplot(aes(x = week, y = mean_log_cd4, color = factor(treatment))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) + 
  labs(x = "Week", y = "Log(CD4)")

## Imbalances in the data 

aids %>%
  group_by(id) %>% 
  count() %>% 
  summary()

aids %>%
  group_by(treatment) %>% 
  distinct(id) %>% 
  count()

# Mistimed events
aids %>% 
  filter(treatment == 2) %>% 
  ggplot(aes(x = week)) +
  geom_histogram(color = "black", fill = "gray") +
  labs(title = "Week of Treatment 2 group")

## Outliers 
aids %>% 
  ggplot(aes(y = log_cd4)) +
  geom_boxplot()

# outliers values
tibble(outliers = boxplot.stats(aids$log_cd4)$out )

# number of outliers 
observerd = tibble(r = boxplot.stats(aids$log_cd4)$out) %>% 
  nrow()

# expected number of outliers are
expected = nrow(aids) * .05

aids %>% 
  ggplot(aes(x = log_cd4)) +
  geom_histogram(bins = 15, color = "black", fill = "gray") +
  labs(title = "Histogram of Log_cd4")

# Standardized of log_cd4 with mean 0 and sd 1 
log_cd4 = aids$log_cd4

scaled = scale(log_cd4)
colMeans(scaled)
apply(scaled, 2, sd)

tibble(r = scaled) %>% 
  ggplot(aes(x = r)) +
  geom_histogram( color = "black", fill = "gray") +
  labs(title = "Standardized Histogram of Log_cd4")


# GLM Modeling ----------------------------------------------------------------

aids = aids %>% 
  mutate(occasion = ceiling(week),
         gender = factor(gender, level = c("male", "female")),
         treatment = as.factor(treatment))


model1 = lme(log_cd4 ~ occasion + I(occasion^2) + 
               treatment + treatment:occasion + age + gender, data = aids,
             random = ~ occasion + I(occasion^2)| id,
             method = "REML")

summary(model1)

model1 = lme(log_cd4 ~ occasion + I(occasion^2) + 
               treatment:occasion + age + gender, data = aids,
             random = ~ occasion + I(occasion^2)| id,
             method = "REML")

summary(model1)

# Including or not including gender
model_no_gender = lme(log_cd4 ~ occasion + I(occasion^2) + 
                        treatment:occasion + treatment:I(occasion^2) +
                        age, data = aids,
                      random = ~ occasion + I(occasion^2)| id,
                      method = "ML")

model = lme(log_cd4 ~ occasion + I(occasion^2) + 
              treatment:occasion + treatment:I(occasion^2) + age + gender, data = aids,
            random = ~ occasion + I(occasion^2)| id,
            method = "ML")

anova(model_no_gender, model)

# Including or not including interactions
model_interactions = lme(log_cd4 ~  occasion + I(occasion^2) + 
                           treatment:occasion + treatment:I(occasion^2) + 
                           age, data = aids,
                         random = ~ occasion + I(occasion^2)| id,
                         method = "ML")

model_without_interaction = lme(log_cd4 ~  occasion + I(occasion^2) + 
                                  age, data = aids,
                                random = ~ occasion + I(occasion^2)| id,
                                method = "ML")

anova(model_interactions, model_without_interaction)

# Linear and Quadratic
model_linear = lme(log_cd4 ~ occasion  + 
                     treatment:occasion +  age, data = aids,
                   random = ~ occasion | id,
                   method = "ML")

model_quadratic = lme(log_cd4 ~ occasion + I(occasion^2) + 
                        treatment:occasion + treatment:I(occasion^2)
                      +  age, data = aids,
                      random = ~ occasion + I(occasion^2)| id,
                      method = "ML")

summary(model_quadratic)

anova(model_linear, model_quadratic)

# Spline
aids <- aids %>% 
  mutate(knot_term = if_else(occasion > 20, occasion, 0)) %>% 
  relocate(knot_term, .after = occasion)

ctrl <- lmeControl(opt = 'optim')
model_spline = lme(log_cd4 ~ occasion + 
                     treatment:occasion + age + knot_term, data = aids,
                   random = ~ occasion + knot_term| id,
                   method = "ML",
                   control = ctrl)
summary(model_spline)

anova(model_linear, model_spline)

# Quadratic Before Knot Term, Linear After Knot Term Spline
aids <- aids %>% 
  mutate(knot_term1 = if_else(treatment == 4, occasion^2, 0)) %>% 
  relocate(knot_term1, .after = occasion)

ctrl <- lmeControl(opt = 'optim')
model_spline2 = lme(log_cd4 ~ occasion + 
                      treatment:occasion + age + knot_term1, data = aids,
                    random = ~ occasion | id,
                    method = "ML",
                    control = ctrl)
summary(model_spline2)

ctrl <- lmeControl(opt = 'optim')
model_spline1 = lme(log_cd4 ~ occasion + 
                      treatment:occasion + age + knot_term1, data = aids,
                    random = ~ occasion + knot_term1| id,
                    method = "ML",
                    control = ctrl)
summary(model_spline1)

anova(model_spline, model_spline1)
anova(model_linear, model_spline1)
anova(model_quadratic, model_spline1)

# Quadratic Spline Model
ctrl <- lmeControl(opt = 'optim')
model_quad_splines = lme(log_cd4 ~ occasion + I(occasion^2) + 
                           treatment:occasion + age + knot_term1, data = aids,
                         random = ~ occasion + I(occasion^2) + knot_term1| id,
                         method = "ML",
                         control = ctrl)
summary(model_quad_splines)
anova(model_quadratic, model_quad_splines)



# Residual Analysis -------------------------------------------------------

res_prog = residuals(model_quadratic, level = 0)


sigma_i = extract.lme.cov(model_quadratic, aids)

#lower triangular matrix
L_i = t(chol(sigma_i))

#transformed residuals
res_prog_trans = solve(L_i) %*% res_prog


tibble(r = res_prog_trans) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(aes(y = stat(density)), bins = 14 ) +
  geom_function(fun = dnorm, color = "blue")

tibble(r = res_prog) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(aes(y = stat(density)), bins = 14 ) +
  geom_function(fun = dnorm, color = "blue")

tibble(r = res_prog_trans) %>% 
  ggplot(aes(sample = r)) + 
  geom_qq_line(color = "blue") +
  geom_qq()

# Transformed Predicted values vs. Transformed Residuals 
mu_hat = fitted(model_quadratic, level = 0)
mu_hat_transformed = solve(L_i) %*% mu_hat


tibble(x = mu_hat_transformed, y = res_prog_trans) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(shape = 1) + 
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Transformed Predicted Value", y = "Transformed Residual")

# Transformed Predicted values vs. Absolute Transformed Residuals  
res_prog_trans_abs = abs(res_prog_trans)

tibble(x = mu_hat_transformed, y = res_prog_trans_abs) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_hline(yintercept = 0.8, linetype = "dashed") + 
  geom_point(shape = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Transformed Predicted Value", y = "Absolute Transformed Residual")

# Mahalanobis Data  
mahalanobis_data = tibble(id = aids$id, r_star = res_prog_trans) %>% group_by(id) %>% 
  nest()

mahalanobis_data = mahalanobis_data %>% 
  mutate(df = map_dbl(data, ~nrow(.x)))

mahalanobis_dist = function(x){
  x = as.matrix(x)
  t(x) %*% x
}

mahalanobis_data = mahalanobis_data %>% 
  mutate(d = map_dbl(data, ~mahalanobis_dist(.x)))

mahalanobis_data = mahalanobis_data %>% 
  mutate(p_value = pchisq(d, df, lower.tail = FALSE))


mahalanobis_data_p = mahalanobis_data %>% 
  arrange(p_value)

mahalanobis_data_p %>% filter(p_value <=.05)

#expected outliers are
expected = 5036 *.05
expected

# Semi-Variogram  
aids = aids %>% 
  mutate(occasion_sqr = occasion ^ 2)

Variogram(model_quadratic,
          data = aids,
          form = ~ occasion + occasion_sqr| id,
          resType = "normalized") %>% 
  as_tibble() %>% 
  ggplot(aes(x = dist, y = variog)) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(shape = 1) +
  geom_smooth(method = "loess", se = FALSE, span = .2)


# GEE and GLME ------------------------------------------------------------

model = glmer(counts2 ~ occasion + I(occasion ^ 2) + treatment:occasion + treatment:I(occasion^2) + age + (1 | id), 
              data = aids,
              family = poisson,
              control = glmerControl(tol = 1e-12),
              nAGQ = 0,
              na.action = na.omit)
summary(model)

model1 = glmer(counts2 ~ occasion + I(occasion ^ 2) + treatment:occasion + treatment:I(occasion^2) + age + (1 + occasion | id), 
               data = aids,
               family = poisson,
               control = glmerControl(tol = 1e-12),
               nAGQ = 0,
               na.action = na.omit)
summary(model1)

anova(model, model1)

model_quadratic2 = lme(log_cd4 ~ occasion + I(occasion^2) + 
                         treatment:occasion + treatment:I(occasion^2)
                       +  age, data = aids,
                       random = ~ occasion + I(occasion^2) | id,
                       method = "ML")
summary(model_quadratic2)

# GLME AIC = 52553.6   
# LME AIC = 11888.81 

