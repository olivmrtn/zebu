make_trial_dataset = function() {

  set.seed(63) # Set seed for reproducibility

  # Simulation parameters
  n <- 100 # Number of patients
  pr <- 0.8 # Proportion of resistant patients
  pe <- 0.2 # Placebo effect on biomarker levels
  de <- 0.3 # Drug effect on biomarker levels
  sd <- 0.05 # Patient and drug effect individual differences (noise)
  sickbiom <- 0.3 # Mean of biomarker level for sick patients
  healthybiom <- 0.7 # Value above which patients are considered healthy

  # Randomly allocated patients to control or test group
  drug <- sample(c("placebo", "drug"), size = n, prob = c(0.5, 0.5), replace = TRUE)

  # Randomly allocated resistance status to patients
  resistance <- sample(c("sensitive", "resistant"), size = n, prob = c(1-pr, pr), replace = TRUE)

  # Randomly allocated unhealthy biomarker levels to patients: pretreatment biomarker levels
  # Simulated with a normal distribution centred around 0.3 with a standard deviation of 0.1
  prebiom <- stats::rnorm(n, mean = sickbiom, sd = sd)
  prebiom <- pmax(pmin(prebiom, 1), 0)

  # Simulate treatment of patients according to their resistance, drug intake and pretreatment biomarker levels
  postbiom <- prebiom + pe + de * (resistance == "sensitive") * (drug == "drug") + stats::rnorm(n, mean = 0, sd = sd)
  postbiom <- pmax(pmin(postbiom, 1), 0)

  # Summarize results in data.frame and save
  trial <- data.frame(drug, resistance, prebiom, postbiom)
  devtools::use_data(trial)
}
