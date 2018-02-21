import os
import csv
import math

wells = []
output = [['well', 'Pwh', 'Ppr', 'Tpr', 'H', 'P_ideal', 'Z', 'P_real', 'Pbh']]
outFile = "results.csv"

# Get the type of program: Real or Ideal
RoI = input("Is this an ideal case? press Y/n")

# Get the wells and their parameters from wells.csv
with open('wells.csv') as File:
    reader = csv.DictReader(File)
    for row in reader:
        wells.append(row)

# Check for the type of program to run
if RoI == "Y" or RoI == "y":
    # For Ideal case the z = 1
    z = 1
elif RoI == "N" or RoI == "n":
    if os.path.exists(outFile):
        z = []
        with open(outFile, 'rU') as File:
            reader = csv.DictReader(File)
            for row in reader:
                z.append(row)
    else:
        print('The file results.csv doesn\'t exist')
        exit()
    
else:
    exit()
# Get the a 
def a(w, Z=1):
    Qgs = float(w["gor"]) * float(w["Qo"])
    return ( (15.36 * float(w["Qs"]) * float(w["ys"])) + (84.93 * ( (float(w["yw"]) * float(w["Qw"])) + (float(w["Qo"]) * float(w["yo"])) )) + ( 0.0188 * float(w["yg"]) * Qgs) ) / ( (float(w["T"])+459.67) * Qgs * Z)

# Get the b
def b(w, Z=1):
    Qgs = float(w["gor"]) * float(w["Qo"])
    return ( (0.025 * float(w["Qs"])) + ( 1.38 * (float(w["Qw"]) + float(w["Qo"])) ) ) / ( (float(w["T"])+459.67) * Qgs * Z )

# Get the c
def c(w, Z=1):
    Qgs = float(w["gor"]) * float(w["Qo"])
    A = (3.141592 * (float(w["P_i"])**2)) / 4
    return ( 0.00678 * (float(w["T"])+459.67) * Qgs * Z ) / A

# Get the d
def d(w, Z=1):
    A = (3.141592 * (float(w["P_i"])**2)) / 4
    return (0.00936 * (float(w["Qw"]) + float(w["Qo"]))) / A

# Get the e
def e(w, Z=1):
    ra = ( 2 * float(w["ep"]) ) / float(w["P_i"])
    f = ( 1 / ( 1.74 - ( 2 * math.log10(ra) ) ) )**2
    return f / ( 2 * 32.2 * float(w["P_i"]) )

            
# Newton Raphson's method function accepts the function, its derivative, starting point for x and epsalum
def Newton(f, dfdx, x, eps):
    f_value = f(x)
    iteration_counter = 0
    while abs(f_value) > eps and iteration_counter < 100:
        try:
            x = x - float(f_value)/dfdx(x)
        except ZeroDivisionError:
            print("Error! - derivative zero for x = ", x)
            sys.exit(1)     # Abort with error

        f_value = f(x)
        iteration_counter += 1
                                      
    # Here, either a solution is found, or too many iterations
    if abs(f_value) > eps:
        iteration_counter = -1
    return x

for well in wells:
    if type(z) is list:
        for row in z:
            if well["well"] == row["well"]:
                Z = float(row["Z"])
    else:
        Z = 1
    # Function to find P1
    def f1(x):
        return (144*b(well, Z)* x) + math.log(x) - ( (144*b(well, Z)*float(well["pwh"])) + math.log(float(well["pwh"])) + ((a(well, Z) + (a(well, Z)*(d(well, Z)**2)*e(well, Z))) * float(well["H"])))

    # The differential of the function f1
    def df1dx(x):
        return (144*b(well, Z)) + (1/x)

    # The function to find P2
    def f2(x):
        return (995328*b(well, Z)*(x**3)) + (10368*(x**2)) - ( (a(well, Z)*(c(well, Z)**2)*e(well, Z)*float(well["H"])) + ( (b(well, Z)*(float(well["pwh"])**3)/3) + (float(well["pwh"])**2)/2 ))

    # The differential of the function f2
    def df2dx(x):
        return (3*(995328*b(well, Z))*(x**2)) + (2*(10368)*x)

    # The function to find P3
    def f3(x):
        return (5184*b(well, Z)*(x**2)) + (72*x) - ( (a(well, Z)*c(well, Z)*d(well, Z)*e(well, Z)*float(well["H"])) + ((b(well, Z)*(float(well["H"])**2))/4) + (float(well["pwh"])/2) )
    
    # The differential of the function f3
    def df3dx(x):
        return (2*5184*b(well, Z)*x) + 72

    P1 = Newton(f1, df1dx, x=1, eps=1.0e-6)  # Use newton's method to get p1

    P2 = Newton(f2, df2dx, x=1, eps=1.0e-6)  # Use Newton's method to get P2

    P3 = Newton(f3, df3dx, x=1, eps=1.0e-6)  # Use Newton's method to get P3
    
    P = P1 + P2 + P3

    PG = P - min([P1, P2, P3])
    
    Tpc = 168 + 325*(float(well["yg"])) + 12.5*(float(well["yg"])**2)

    Ppc = 677 + 15.0*(float(well["yg"])) + 37.5*(float(well["yg"])**2)

    Tpr = (float(well["T"])+459.67) / Tpc

    Ppr = P / Ppc

    if type(z) is int:
        output.append([well["well"], well["pwh"], Ppr, Tpr, well["H"], P, '', '', PG])
    else:
        for row in z:
            if row["well"] == well["well"]:
                output.append([row["well"], row["Pwh"], row["Ppr"], row["Tpr"], row["H"], row["P_ideal"], row["Z"], P, row["Pbh"]])
                stout = "Name of Well: {0} \n\n Well Head Pressure Pwh (Psia) = {1} \n Bottom Hole Temperature T = {2} \n Well Depth H (ft) = {3} \n Tubing Internal Diameter ID (inches) = {4} \n Roughness Factor = {5} \n Compressibility Factor for ideal case = {6} \n Compressibility Factor for Real case = {7} \n P1 = {8} \n P2 = {9} \n P3 = {10} \n Calculated bottom hole Pressue for Ideal case = {11} \n Calculated bottom hole Pressure for Real case = {12} \n\n\n".format(row["well"], float(well['pwh']), float(well["T"])+459.67, float(well["H"]), float(well["P_i"]), float(well["ep"]), 1, float(row["Z"]), P1, P2, P3, float(row["P_ideal"]), P)
                print(stout)

with open('results.csv', 'w') as File:
    writer = csv.writer(File)      
    writer.writerows(output)

input("Press Enter to continue...")
