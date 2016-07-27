set.seed(63) # Set seed for reproducibility

# Simulation parameters
n <- 100 # Number of patients
pr <- 0.5 # Proportion of resistant patients
pe <- 0.3 # Placebo effect on biomarker levels
de <- 0.3 # Drug effect on biomarker levels
sd <- 0.05 # Patient and drug effect individual differences
sickbiom <- 0.3 # Mean of biomarker level for sick patients
healthybiom <- 0.7 # Value above which patients are considered healthy

# Randomly allocated patients to control or test group
drug <- sample(c("placebo", "drug"), size = n, prob = c(0.5, 0.5), replace = TRUE)

# Randomly allocated resistance status to patients
resistance <- sample(c("sensitive", "resistant"), size = n, prob = c(1-pr, pr), replace = TRUE)

# Randomly allocated unhealthy biomarker levels to patients: pretreatment biomarker levels
# Simulated with a normal distribution centred around 0.3 with a standard deviation of 0.1
prebiom <- rnorm(n, mean = sickbiom, sd = sd)

# Simulate treatment of patients according to their resistance, drug intake and pretreatment biomarker levels
postbiom <- prebiom + pe + de * (resistance == "sensitive") * (drug == "drug") + rnorm(n, mean = 0, sd = sd)

# Summarize results in data.frame and save
trial <- data.frame(drug, resistance, prebiom, postbiom)
devtools::use_data(trial)

las <- lassie(trial, select = c("drug", "postbiom"), continuous = "postbiom", breaks = c(0, healthybiom, 1))
las <- permtest(las)
plot(las, "local")
plot(las, "obs")

sub <- subgroups(x = las, y = trial, select = "resistance", alpha = 0.01)
sub <- permtest(sub, nb = 1000)
plot(sub, "local")
plot(sub, "obs")

x = table(drug = trial$drug, healthy = trial$postbiom > 0.6)
x
