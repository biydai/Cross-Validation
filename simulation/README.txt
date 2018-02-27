This file lists all simulation scenarios that were considered in the cross-validation paper

# Scenario 1 (plots the lambda selection and MSE ratio)
compare four cross validation methods: 
ungrouped, grouped, linear predictor, deviance residuals
varying magnitude of betas at higher dimension (n = 200, p = 1000, nonzero = 10)

# Scenario 1.2 (plots the lambda selection and MSE ratio)
compare four cross validation methods: 
ungrouped, grouped, linear predictor, deviance residuals
varying magnitude of betas at low dimension (n = 100, p = 10, nonzero = 2)

# Scenario 1.3 (plots the lambda selection and MSE ratio)
compare four cross validation methods: 
ungrouped, grouped, linear predictor, deviance residuals
varying magnitude of betas at low dimension (n = 100, p = 100, nonzero = 10)

# Scenario 1.4 (plots the lambda selection and MSE ratio)
compare four cross validation methods: 
ungrouped, grouped, linear predictor, deviance residuals
varying magnitude of betas at low dimension (n = 100, p = 1000, nonzero = 10)

# Scenario 2.1 (put in tables with lambda(SE), MSE)
compare four cross validation methods
at low dimension (n = 100, p = 20, nonzero = 10)

# Scenario 2.2
compare four cross validation methods
at low dimension (n = 100, p = 100, nonzero = 10)

# Scenario 2.3
compare four cross validation methods
at low dimension (n = 200, p = 1000, nonzero = 10)

# Scenario 2.4 
compare four cross validation methods
at high dimension (n = 500, p = 10,000, nonzero = 20) 

# Scenario 3.1
compare four cross validation methods in terms of stability
varying censored portion

# Scenario 3.2
compare four cross validation methods in terms of stability
varying number of folds

# Scenario 4.1
compare four cross validation methods with AIC, BIC, EBIC?, GCV?
varying censored portion

# Scenario 5
Change the penalty to MCP, just to see what has happened
Repeat 2.1,2.2 and 2.3 

# Scenario 6
Change baseline estimation for deviance residuals
exact
KP's method
weighted KP's method
Breslow's method
parametric approach