# This file contains functions that might be helpful for Exam 2 of EEL3111C for the Fall '24 semester
#   PRESS "CTRL+R" TO RUN!
#   Contents
#       - rlc_response()
#           > Inputs: a circuit's resistance, capacitance, and inducatance and if it is parallel/series and if looking for natural/step response
#           > Outputs: whether overdamped, underdamped, critically damped, roots of char eq, neper frequency, resonant frequency, all associated equations
#       - Notes for filter-related problems
#       - Miscellaneous notes
#
#   Important:
#       The TI-NSpire runs an older version of python, and therefore does not support the following python elements
#           > Match/Case statements (use if/elif/else instead)
#           > fstrings (alternative example: myStr = "value1 = {}, value2 = {}".format(v1, v2) )
#           > __name__ == __main__ (just put code to be automatically run at the very end of the .py file outside any functions)

# <Settings>
verboseOutput = False
# </Settings>

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
    

# Function for natural/step response for parallel/series RLC circuits
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
        charEq = "s^2 + {:.2e}*s + {:.2e} = 0".format(1/(r*c), 1/(l*c))

        # Set up equations v(t) [v_t], v(0+) [v_t_0], dv(0+)/dt [deriv]
        if state == "overdamped":
            eq1 = "v(t) = A1*e^({:.2e}*t) + A2*e^({:.2e}*t), t>=0".format(s1, s2)
            eq2 = "v(0+) = A1 + A2 = V0"
            eq3 = "dv(0+)/dt = {:.2e}*A1 + {:.2e}*A2 = -(V0 + R*I0)/(R*C)".format(s1, s2)
        elif state == "underdamped":
            eq1 = "v(t) = B1*e^({:.2e}*t)*cos({:.2e}*t) + B2*e^({:.2e}*t)*sin({:.2e}*t), t>=0".format(-a,w,-a,w)
            eq2 = "v(0+) = B1 = V0"
            eq3 = "dv(0+)/dt = {:.2e}*B1 + {:.2e}*B2 = -(V0 + R*I0)/(R*C)".format(-a,w)
        elif state == "critically damped":
            eq1 = "v(t) = D1*t*e^({:.2e}*t) + D2*e^({:.2e}*t), t>=0".format(-a,-a)
            eq2 = "v(0+) = D2 = V0"
            eq3 = "dv(0+)/dt = D1 - {:.2e}*D2 = -(V0 + R*I0)/(R*C)".format(a)



    if ( (orientation == "p") and (typeResponse == "s") ): #parallel step response

        # Calculate neper freq (a) and resonant freq (w)
        a = 1/(2*r*c)
        w = math.sqrt(1/(l*c))

        # Determine damped state
        state = getDampedState(a,w)

        # Calculate roots of characteristic equation
        s1,s2 = calculateRoots(a,w)

        # Set up characteristic equation
        charEq = "s^2 + {:.2e}*s + {:.2e} = 0".format(1/(r*c), 1/(l*c))

        # Set up equations iL(t) [iL_t], iL(0+) [iL_t_0], diL(0+)/dt [deriv] 
        if state == "overdamped":
            eq1 = "iL(t) = I_f + A1*e^({:.2e}*t) + A2*e^({:.2e}*t), t>=0".format(s1, s2)
            eq2 = "iL(0+) = I_f + A1 + A2 = I_0"
            eq3 = "diL(0+)/dt = {:.2e}*A1 + {:.2e}*A2 = (V_0/L)".format(s1, s2)
        elif state == "underdamped":
            eq1 = "iL(t) = I_f + B1*e^({:.2e}*t)*cos({:.2e}*t) + B2*e^({:.2e}*t)*sin({:.2e}*t), t>=0".format(-a, w, -a, w)
            eq2 = "iL(0+) = I_f + B1 = I_0"
            eq3 = "diL(0+)/dt = {:.2e}*B1 + {:.2e}*B2 = (V_0/L)".format(-a, w)
        elif state == "critically damped":
            eq1 = "iL(t) = I_f + D1*t*e^({:.2e}*t) + D2*e^({:.2e}*t), t>=0".format(-a,-a)
            eq2 = "iL(0+) = I_f + D2 = I_0"
            eq3 = "diL(0+)/dt = D1 - {:.2e}*D2 = (V_0/L)".format(a)


    if ( (orientation == "s") and (typeResponse == "n") ): #series natural response

        # Calculate neper freq (a) and resonant freq (w)
        a = r/(2*l)
        w = 1/math.sqrt(l*c)

        # Determine damped state
        state = getDampedState(a,w)

        # Calculate roots of characteristic equation
        s1,s2 = calculateRoots(a,w)

        # Set up characteristic equation
        charEq = "s^2 + {:.2e}*s + {:.2e} = 0".format(r/l, 1/(l*c))

        # Set up equations iL(t) [iL_t], iL(0+) [iL_t_0], diL(0+)/dt [deriv]
        if state == "overdamped":
            eq1 = "iL(t) = A1*e^({:.2e}*t) + A2*e^({:.2e}*t), t>=0".format(s1, s2)
            eq2 = "iL(0+) = A1 + A2 = I_0"
            eq3 = "diL(0+)/dt = {:.2e}*A1 + {:.2e}*A2 = -(1/L)*(R*I_0 + V_0)".format(s1,s2)
        elif state == "underdamped":
            eq1 = "iL(t) = B1*e^({:.2e}*t)*cos({:.2e}*t) + B2*e^({:.2e}*t)*sin({:.2e}*t), t>=0".format(-a,w,-a,w)
            eq2 = "iL(0+) = B1 = I_0"
            eq3 = "diL(0+)/dt = {:.2e}*B1 + {:.2e}*B2 = -(1/L)*(R*I_0 + V_0)".format(-a,w)
        elif state == "critically damped":
            eq1 = "iL(t) = D1*t*e^({:.2e}*t) + D2*e^({:.2e}*t), t>=0".format(-a,-a)
            eq2 = "iL(0+) = D2 = I_0"
            eq3 = "diL(0+)/dt = D1 - {:.2e}*D2 = -(1/L)*(R*I_0 + V_0)".format(a)
                

    
    if ( (orientation == "s") and (typeResponse == "s") ): #series step response

        # Calculate neper freq (a) and resonant freq (w)
        a = r/(2*l)
        w = 1/math.sqrt(l*c)

        # Determine damped state
        state = getDampedState(a,w)

        # Calculate roots of characteristic equation
        s1,s2 = calculateRoots(a,w)

        # Set up characteristic equation
        charEq = "s^2 + {:.2e}*s + {:.2e} = V/{:.2e}".format(r/l, 1/(l*c), l*c)

        # Set up equations vC(t) [vC_t], vC(0+) [vC_t_0], dvC(0+)/dt [deriv]
        if state == "overdamped":
            eq1 = "vC(t) = V_f + A1*e^({:.2e}*t) + A2*e^({:.2e}*t), t>=0".format(s1, s2)
            eq2 = "vC(0+) = V_f + A1 + A2 = V_0"
            eq3 = "dvC(0+)/dt = {:.2e}*A1 + {:.2e}*A2 = (I_0/C)".format(s1, s2)
        elif state == "underdamped":
            eq1 = "vC(t) = V_f + B1*e^({:.2e}*t)*cos({:.2e}*t) + B2*e^({:.2e}*t)*sin({:.2e}*t), t>=0".format(-a, w, -a, w)
            eq2 = "vC(0+) = V_f + B1 = V_0"
            eq3 = "dvC(0+)/dt = {:.2e}*B1 + {:.2e}*B2 = (I_0/C)".format(-a,w)
        elif state == "critically damped":
            eq1 = "vC(t) = V_f + D1*t*e^({:.2e}*t) + D2*e^({:.2e}*t), t>=0".format(-a, -a)
            eq2 = "vC(0+) = V_f + D2 = V_0"
            eq3 = "dvC(0+)/dt = D1 - {:.2e}*D2 = (I_0/C)".format(a)

    # Print results
    if (orientation == "p"):
        orien_str = "Parallel"
    else:
        orien_str = "Series"

    if (typeResponse == "p"):
        type_str = "Natural"
    else:
        type_str = "Step"

    print("\n\nOutputs:")
    print("\t>Circuit is {}".format(state))
    print("\t>Neper freq = {:.2e}".format(a))
    print("\t>Resonant freq = {:.2e}".format(w))

    if (state == "underdamped"):
        w_d = math.sqrt(w**2 - a**2)
        print("\tDamped freq = {:.2e}".format(w_d))

    print("\t>s1 = {}".format(s1))
    print("\t>s2 = {}".format(s2))

    # The equations for each config only print if global variable verboseOutput is true

    if (verboseOutput == True):
        print("\t>Char Eq: {}".format(charEq))
        print("\t>{}".format(eq1))
        print("\t>{}".format(eq2))
        print("\t>{}".format(eq3))

    return

# Function for use with loaded filter problems
def loaded_filter():
    print("What filter config do you have?")
    print("\t1) low-pass RL")
    print("\t2) high-pass RL")
    print("\t3) low-pass RC")
    print("\t4) high-pass RC")

    option = input("Config:" )

    print("\n")

    if (option == '1'): #low-pass RL
        print("Solve for load R_eq, then use")
        print("\tcutoff freq w = R/L")
    elif (option == '2'): #high-pass RL
        print("Equation to use:")
        print("\tR_L = (w*L * R) / (R - w*L)")
        print()
        print("\tR_L: load res resistance")
        print("\tw: cutoff freq (rads/s)")
    elif (option == '3'):
        print("Equation to use:")
        print("\tR_L = (R) / (w*C*R - 1)")
        print()
        print("\tR_L: load res resistance")
        print("\tw: cutoff freq (rads/s)")
    elif (option == '4'): #high-pass RC
        print("Solve for load R_eq, then use")
        print("\tcutoff freq w = 1/(R*C)")

    return

# Function for all things related to RC/RL/RLC passive filter design
def filter_design():
    print("\nChoose an action:")
    print("\t1) Print cap/ind impedance formulas")
    print("\t2) How to identify high/low pass?")
    print("\t3) Print Vout and |Vout| formulas")
    print("\t4) Loaded filter problem")
    print("\t5) Determining the filter gain")
    print("\t6) Creating a Bode plot")

    option = input("Selection: ")

    if (option == '1'): #printing capacitor and inductor impedance formulas
        print("\nImpedance Formulas (z)")
        print("\tZ_c = 1 / (w * c * j)")
        print("\tZ_l = w * L * j")

    elif (option == '2'): #printing tips for identifying a high or low-pass RC/RL filter
        print("\nIdentifing high/low pass filters")
        print("\tLow pass:")
        print("\t\tInductor in series w/ src")
        print("\t\t(Low freqs see a short)")
        print("\t\tCap to gnd")
        print("\tHigh pass:")
        print("\t\tInductor to gnd")
        print("\t\tCap in series w/ src")
        print("\t\t(High freqs see a short)")

    elif (option == '3'): #printing formulas for Vout and |Vout| for high or low-pass RC/RL filters
        print("\nPrinting Vout & |Vout| formulas")
        print("\tFor low-pass RL filter")
        print("\t\tVout = (R * Vin) / (R + w*L)")
        print("\t\t|Vout| = (R*|Vin|) / sqrt[(R)^2 + (w*L)^2)]")

        print("\tFor high-pass RL filter")
        print("\t\tVout = (Vin * w*L) / (R + w*L))")
        print("\t\t|Vout| = (w*L*|Vin|) / (sqrt[R^2 + (w*L)^2])")

        print("\tFor low-pass RC filter")
        print("\t\tVout = (1 * Vin) / (1 + w*R*C)")
        print("\t\t|Vout| = (1*|Vin|) / sqrt[1+(w*R*C)^2]")

        print("\tFor high-pass RC filter")
        print("\t\tVout = (Vin * w*R*C) / (1 + w*R*C)")
        print("\t\t|Vout| = (w*R*C*|Vin|) / sqrt[1+(w*R*C)^2]")

    elif (option == '4'): #loaded filter problems
        loaded_filter()

    elif (option == '5'): #determining filter gain
        print("Determining filter gain")
        
        print("\tLow-pass RL:")
        print("\t\t|H(w)| = w*L / sqrt[R^2 + (w*L)^2]")
        
        print("\tHigh-pass RL:")
        print("\t\t|H(w)| = R / sqrt[R^2 + (w*L)^2]")
        
        print("\tLow-pass RC:")
        print("\t\t|H(w)| = 1 / sqrt[1 + (w*R*C)^2]")
        
        print("\tHigh-pass RC:")
        print("\t\t|H(w)| = (w*R*C) / sqrt[1 + (w*R*C)^2]")

    elif (option == '6'): #creating a Bode plot
        print("Creating a Bode plot")
        print("\tEq: 20*log(10)(H(jw))")
        print()
        print("\tCutoff freq @ -3dB")
        print("\tGain of 1 @ 0dB")

    return
    
# Function for miscellaneous tips
def misc_tips():
    print("Topic in question?")
    print("\t1) Unit step input")

    option = input("Topic: ")
    print()

    if (option == 1):
        print("Unit step input")
        print("\tEq: 1-e^(-t/rho)")
        print()
        print("For t=rho,")
        print("1-e^(-1) = 63.2%")

    return

# *****Code to execute when this script is run*****
while(1):
    print("\n\n")
    print("Select a function to use:")
    print("\t1) rlc_response")
    print("\t2) Filter design questions")
    print("\t3) Misc tips")
    print("\tz) Exit program")

    option = input("Function: ")

    if option == '1':
        rlc_response()
    elif option == '2':
        filter_design()
    elif option == '3':
        misc_tips()
    elif option == 'z':
        break
