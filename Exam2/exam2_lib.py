# This file contains functions that might be helpful for Exam 2 of EEL3111C for the Fall '24 semester
#   Contents
#       - rlc_response()
#           > Inputs: a circuit's resistance, capacitance, and inducatance and if it is parallel/series and if looking for natural/step response
#           > Outputs: whether overdamped, underdamped, critically damped, roots of char eq, neper frequency, resonant frequency, all associated equations

import math
import cmath

# Function to determine damped state (used only by rlc_response())
def getDampedState(a,w): #this returns "overdamped, underdamped", or "critically damped" depending on the values of a and w
    if (a > w):
        return "overdamped"
    elif (a < w):
        return "underdamped"
    else:
        return "critically damped"
    

# Function to calculate roots of characteristic equation (used only by rlc_response())
def calculateRoots(a,w):
    try:    # get root s1
        s1 = -a + math.sqrt(a**2 - w**2)
    except ValueError:
        s1 = -a + cmath.sqrt(a**2 - w**2)

    try: #get root s2
        s2 = -a - math.sqrt(a**2 - w**2)
    except ValueError:
        s2 = -a - cmath.sqrt(a**2 - w**2)

    return s1,s2
    

# Function  for natural/step response for parallel/series RLC circuits
def rlc_response():
    # Get parallel or series
    orientation = input("Parallel or Series RLC? (p/s): ")

    # Get natural or step response
    typeResponse = input("Natural or Step Response? (n/s): ")
    
    # Get values
    r = input("Input equivalent resistance (in ohms): ")
    c = input("Input capacitance (in farads): ")
    l = input("Input inductance (in henries): ")

    # Convert values to floats
    r = float(r)
    c = float(c)
    l = float(l)


    # Continue calculations based on circuit orientation and type of response
    if ( (orientation == "p") and (typeResponse == "n") ): #parallel natural response

        # Calculate neper freq (a) and resonant freq (w)
        a = 1/(2*r*c)
        w = 1/math.sqrt(l*c)

        # Determine damped state
        state = getDampedState(a,w)

        # Calculate roots of characteristic equation
        s1,s2 = calculateRoots(a,w)

        # Set up characteristic equation
        charEq = f"s^2 + {1/(r*c):.2e}*s + {1/(l*c):.2e} = 0"

        # Set up equations v(t) [v_t], v(0+) [v_t_0], dv(0+)/dt [deriv]
        match state:
            case "overdamped":
                eq1 = f"v(t) = A1*e^({s1:.2e}*t) + A2*e^({s2:.2e}*t), t>=0"
                eq2 = f"v(0+) = A1 + A2 = V0"
                eq3 = f"dv(0+)/dt = {s1:.2e}*A1 + {s2:.2e}*A2 = -(V0 + R*I0)/(R*C)"
            case "underdamped":
                eq1 = f"v(t) = B1*e^({-1*a:.2e}*t)*cos({w}*t) + B2*e^({-1*a:.2e}*t)*sin({w:.2e}*t), t>=0"
                eq2 = f"v(0+) = B1 = V0"
                eq3 = f"dv(0+)/dt = {-1*a:.2e}*B1 + {w:.2e}*B2 = -(V0 + R*I0)/(R*C)"
            case "critically damped":
                eq1 = f"v(t) = D1*t*e^({-1*a:.2e}*t) + D2*e^({-1*a:.2e}*t), t>=0"
                eq2 = f"v(0+) = D2 = V0"
                eq3 = f"dv(0+)/dt = D1 - {a:.2e}*D2 = -(V0 + R*I0)/(R*C)"


    if ( (orientation == "p") and (typeResponse == "s") ): #parallel step response

        # Calculate neper freq (a) and resonant freq (w)
        a = 1/(2*r*c)
        w = math.sqrt(1/(l*c))

        # Determine damped state
        state = getDampedState(a,w)

        # Calculate roots of characteristic equation
        s1,s2 = calculateRoots(a,w)

        # Set up characteristic equation
        charEq = f"s^2 + {1/(r*c):.2e}*s + {1/(l*c):.2e} = 1/{l*c:.2e}"

        # Set up equations iL(t) [iL_t], iL(0+) [iL_t_0], diL(0+)/dt [deriv] 
        match state:
            case "overdamped":
                eq1 = f"iL(t) = I_f + A1*e^({s1:.2e}*t) + A2*e^({s2:.2e}*t), t>=0"
                eq2 = f"iL(0+) = I_f + A1 + A2 = I_0"
                eq3 = f"diL(0+)/dt = {s1:.2e}*A1 + {s2:.2e}*A2 = (V_0/L)"
            case "underdamped":
                eq1 = f"iL(t) = I_f + B1*e^({-1*a:.2e}*t)*cos({w:.2e}*t) + B2*e^({-1*a:.2e}*t)*sin({w:.2e}*t), t>=0"
                eq2 = f"iL(0+) = I_f + B1 = I_0"
                eq3 = f"diL(0+)/dt = {-1*a:.2e}*B1 + {w:.2e}*B2 = (V_0/L)"
            case "critically damped":
                eq1 = f"iL(t) = I_f + D1*t*e^({-1*a:.2e}*t) + D2*e^({-1*a:.2e}*t), t>=0"
                eq2 = f"iL(0+) = I_f + D2 = I_0"
                eq3 = f"diL(0+)/dt = D1 - {a:.2e}*D2 = (V_0/L)"


    if ( (orientation == "s") and (typeResponse == "n") ): #series natural response

        # Calculate neper freq (a) and resonant freq (w)
        a = r/(2*l)
        w = 1/math.sqrt(l*c)

        # Determine damped state
        state = getDampedState(a,w)

        # Calculate roots of characteristic equation
        s1,s2 = calculateRoots(a,w)

        # Set up characteristic equation
        charEq = f"s^2 + {r/l:.2e}*s + {1/(l*c):.2e} = 0"

        # Set up equations iL(t) [iL_t], iL(0+) [iL_t_0], diL(0+)/dt [deriv] 
        match state:
            case "overdamped":
                eq1 = f"iL(t) = A1*e^({s1:.2e}*t) + A2*e^({s2:.2e}*t), t>=0"
                eq2 = f"iL(0+) = A1 + A2 = I_0"
                eq3 = f"diL(0+)/dt = {s1:.2e}*A1 + {s2:.2e}*A2 = -(1/L)*(R*I_0 + V_0)"
            case "underdamped":
                eq1 = f"iL(t) = B1*e^({-1*a:.2e}*t)*cos({w:.2e}*t) + B2*e^({-1*a:.2e}*t)*sin({w:.2e}*t), t>=0"
                eq2 = f"iL(0+) = B1 = I_0"
                eq3 = f"diL(0+)/dt = {-1*a:.2e}*B1 + {w:.2e}*B2 = -(1/L)*(R*I_0 + V_0)"
            case "critically damped":
                eq1 = f"iL(t) = D1*t*e^({-1*a:.2e}*t) + D2*e^({-1*a:.2e}*t), t>=0"
                eq2 = f"iL(0+) = D2 = I_0"
                eq3 = f"diL(0+)/dt = D1 - {a:.2e}*D2 = -(1/L)*(R*I_0 + V_0)"

    
    if ( (orientation == "s") and (typeResponse == "s") ): #series step response

        # Calculate neper freq (a) and resonant freq (w)
        a = 1/(2*r*c)
        w = 1/math.sqrt(l*c)

        # Determine damped state
        state = getDampedState(a,w)

        # Calculate roots of characteristic equation
        s1,s2 = calculateRoots(a,w)

        # Set up characteristic equation
        charEq = f"s^2 + {r/l:.2e}*s + {1/(l*c):.2e} = V/{l*c:.2e}"

        # Set up equations vC(t) [vC_t], vC(0+) [vC_t_0], dvC(0+)/dt [deriv]
        match state:
            case "overdamped":
                eq1 = f"vC(t) = V_f + A1*e^({s1:.2e}*t) + A2*e^({s2:.2e}*t), t>=0"
                eq2 = f"vC(0+) = V_f + A1 + A2 = V_0"
                eq3 = f"dvC(0+)/dt = {s1:.2e}*A1 + {s2:.2e}*A2 = (I_0/C)"
            case "underdamped":
                eq1 = f"vC(t) = V_f + B1*e^({-1*a:.2e}*t)*cos({w:.2e}*t) + B2*e^({-1*a:.2e}*t)*sin({w:.2e}*t), t>=0"
                eq2 = f"vC(0+) = V_f + B1 = V_0"
                eq3 = f"dvC(0+)/dt = {-1*a:.2e}*B1 + {w:.2e}*B2 = (I_0/C)"
                pass
            case "critically damped":
                eq1 = f"vC(t) = V_f + D1*t*e^({-1*a:.2e}*t) + D2*e^({-1*a:.2e}*t), t>=0"
                eq2 = f"vC(0+) = V_f + D2 = V_0"
                eq3 = f"dvC(0+)/dt = D1 - {a:.2e}*D2 = (I_0/C)"
                pass

    # Print results
    if (orientation == "p"):
        orien_str = "Parallel"
    else:
        orien_str = "Series"

    if (typeResponse == "p"):
        type_str = "Natural"
    else:
        type_str = "Step"

    print("Outputs:")
    print(f"\tCircuit is {state}")
    print(f"\tneper freq = {a:.2e}")
    print(f"\tresonant freq = {w:.2e}")
    print(f"\ts1 = {s1:.2e}")
    print(f"\ts2 = {s2:.2e}")
    print(f"\tChar Eq: {charEq}")
    print(f"\t{eq1}")
    print(f"\t{eq2}")
    print(f"\t{eq3}")

    return


if __name__ == "__main__":
    rlc_response()
