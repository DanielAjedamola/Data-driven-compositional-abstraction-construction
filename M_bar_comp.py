# """
# Approximating the interconnection map M for abstraction construction for 6 tank example
# Created by Daniel A (Oct. '24)
# """

import numpy as np
import itertools
from sklearn.linear_model import LassoCV, Lasso, MultiTaskLassoCV
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt

# obtain center of grid
def generate_grid_center(a, b, cell_size, N):
    point = np.arange(a, b, cell_size)
    centers = []
    center = [x + cell_size/2 for x in point]
    for c in center:
        if (a < c < b):
            centers.append(c)
    point_centers = list(itertools.product(centers, repeat=N))
    return point_centers

N2, N3, N5, N6, N7, N8, N9, N10 = 1, 2, 1, 4, 1, 1, 2, 2
a, b = 0, 100
rho_y, rho_z = 2.5, 5

# Generate grid and collect centers for w_i, i=2,3,5,6
centersy2 = generate_grid_center(a, b, rho_y, N2)
centersy3 = generate_grid_center(a, b, rho_y, N3)
centersy5 = generate_grid_center(a, b, rho_y, N5)
centersy6 = generate_grid_center(a, b, rho_y, N6)
centersy7 = generate_grid_center(a, b, rho_y, N7)
centersy8 = generate_grid_center(a, b, rho_y, N8)
centersy9 = generate_grid_center(a, b, rho_y, N9)
centersy10 = generate_grid_center(a, b, rho_y, N10)
# print(centersy2)



def sq(num):
    return np.sqrt(num)
def sq2(num1,num2):
    return 0.5*(sq(num1)+sq(num2))
def sq4(num1,num2,num3,num4):
    return (sq(num1)+sq(num2)+sq(num3)+sq(num4))/4.0

# interconnection map 
# create a list of (y_i, w_i) tuples
print("Now creating dataset")
Data_sety2 = [(c, sq(c)) for c in centersy2]
Data_sety3 = [(c, sq2(c[0],c[1])) for c in centersy3]
Data_sety5 = [(c, sq(c)) for c in centersy5]
Data_sety6 = [(c, sq4(c[0],c[1],c[2],c[3])) for c in centersy6]
Data_sety7 = [(c, sq(c)) for c in centersy7]
Data_sety8 = [(c, sq(c)) for c in centersy8]
Data_sety9 = [(c, sq2(c[0],c[1])) for c in centersy9]
Data_sety10 = [(c, sq2(c[0],c[1])) for c in centersy10]

# Print the collected centers
print(len(Data_sety2),len(Data_sety3),len(Data_sety5),len(Data_sety6))

# cross validation to tune the most appriopriate lasso parameter lambda
# Assume X is your centers
Xy2, w2 = zip(*Data_sety2)
Xy3, w3 = zip(*Data_sety3)
Xy5, w5 = zip(*Data_sety5)
Xy6, w6 = zip(*Data_sety6)
Xy7, w7 = zip(*Data_sety7)
Xy8, w8 = zip(*Data_sety8)
Xy9, w9 = zip(*Data_sety9)
Xy10, w10 = zip(*Data_sety10)

# Define a range of alpha values to search over
lambdas = np.logspace(-3, 3, 9)  # Example range, adjust as needed

print("Now CVing for the best lambda")
# Define the number of folds for cross-validation
k = 5

# Initialize the cross-validator
kf = KFold(n_splits=k, shuffle=True, random_state=42)

# Initialize the LassoCV model with a range of alphas and perform cross-validation
# lasso_cv = MultiTaskLassoCV(cv=KFold(n_splits=5))
lasso_cvy2 = LassoCV(alphas=lambdas, cv=kf)
lasso_cvy2.fit(Xy2, w2)
lasso_cvy3 = LassoCV(alphas=lambdas, cv=kf)
lasso_cvy3.fit(Xy3, w3)
lasso_cvy5 = LassoCV(alphas=lambdas, cv=kf)
lasso_cvy5.fit(Xy5, w5)
lasso_cvy6 = LassoCV(alphas=lambdas, cv=kf)
lasso_cvy6.fit(Xy6, w6)
lasso_cvy7 = LassoCV(alphas=lambdas, cv=kf)
lasso_cvy7.fit(Xy7, w7)
lasso_cvy8 = LassoCV(alphas=lambdas, cv=kf)
lasso_cvy8.fit(Xy8, w8)
lasso_cvy9 = LassoCV(alphas=lambdas, cv=kf)
lasso_cvy9.fit(Xy9, w9)
lasso_cvy10 = LassoCV(alphas=lambdas, cv=kf)
lasso_cvy10.fit(Xy10, w10)

# Get the best alpha value
best_lambday2 = lasso_cvy2.alpha_
best_lambday3 = lasso_cvy3.alpha_
best_lambday5 = lasso_cvy5.alpha_
best_lambday6 = lasso_cvy6.alpha_
best_lambday7 = lasso_cvy7.alpha_
best_lambday8 = lasso_cvy8.alpha_
best_lambday9 = lasso_cvy9.alpha_
best_lambday10 = lasso_cvy10.alpha_
print(f"Best lambda value: y; {best_lambday2}, {best_lambday3}, {best_lambday5} and {best_lambday6}")

print("Now fitting model with best lambda")
# Initialize and train the Lasso model with the best lambda value
lassoy2 = Lasso(alpha=best_lambday2)
lassoy2.fit(Xy2, w2)
lassoy3 = Lasso(alpha=best_lambday3)
lassoy3.fit(Xy3, w3)
lassoy5 = Lasso(alpha=best_lambday5)
lassoy5.fit(Xy5, w5)
lassoy6 = Lasso(alpha=best_lambday6)
lassoy6.fit(Xy6, w6)
lassoy7 = Lasso(alpha=best_lambday7)
lassoy7.fit(Xy7, w7)
lassoy8 = Lasso(alpha=best_lambday8)
lassoy8.fit(Xy8, w8)
lassoy9 = Lasso(alpha=best_lambday9)
lassoy9.fit(Xy9, w9)
lassoy10 = Lasso(alpha=best_lambday10)
lassoy10.fit(Xy10, w10)

# Get the coefficients (M) from the trained model
My2 = lassoy2.coef_
print("Matrix My2:")

# Print out the coefficients
for i, coef in enumerate(My2):
    print(f"Feature {i+1}: {coef}")

My3 = lassoy3.coef_
print("Matrix My3:")

# Print out the coefficients
for i, coef in enumerate(My3):
    print(f"Feature {i+1}: {coef}")

My5 = lassoy5.coef_
print("Matrix My5:")

# Print out the coefficients
for i, coef in enumerate(My5):
    print(f"Feature {i+1}: {coef}")

My6 = lassoy6.coef_
print("Matrix My6:")

# Print out the coefficients
for i, coef in enumerate(My6):
    print(f"Feature {i+1}: {coef}")

My7 = lassoy7.coef_
print("Matrix My7:")

# Print out the coefficients
for i, coef in enumerate(My7):
    print(f"Feature {i+1}: {coef}")

My8 = lassoy8.coef_
print("Matrix My8:")

# Print out the coefficients
for i, coef in enumerate(My8):
    print(f"Feature {i+1}: {coef}")

My9 = lassoy9.coef_
print("Matrix My9:")

# Print out the coefficients
for i, coef in enumerate(My9):
    print(f"Feature {i+1}: {coef}")

My10 = lassoy10.coef_
print("Matrix My10:")

# Print out the coefficients
for i, coef in enumerate(My10):
    print(f"Feature {i+1}: {coef}")