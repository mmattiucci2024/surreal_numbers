###################################################################################
# File: mysn.py
# Ver: 1.7
# Updated: 10th of November 2025 
# GOAL: This SW provides procedures and fucntions for handling surreal numbers 
#       and their conversion to/from decimals, integers, and fractions. 
#       It also provides procedures for addition, subtraction, sign change, 
#       multiplication, inversion, and division operations between surreal numbers.
#       A surreal number X is a pair of two lists each of which recursively 
#       contains surreal numbers:
#       X = (LEFTX,RIGHTX)
#       LEFTX = [Xl1,Xl2,...,Xln]
#       RIGHTX = [Xr1,Xr2,...,Xrm]
#       where Xl1,Xl2,...Xln and Xr1,Xr2,...,Xrm are surreal numbers such 
#       that it does not exist Xli greater or equal Xrj for any i,j.
###################################################################################
#
# MIT License
#
# Copyright (c) 2025 Marco Mattiucci
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###################################################################################



import math
import random



############################################################################################
# BASIC DEFINITION OF SURREAL ENTITIES:
# 0) the empty set does exist and is: []
# 1) surreal_ZERO() is a couple of two empty sets.
# 2) surreal_ONE() is a couple of a set with ZERO (left) and an empty set (right).
# 3) surreal_NEGATIVE_ONE() is a couple of an empty set (left) and a set with ZERO (right).
# these are axiomatic truths.
############################################################################################

# Definition of surreal ZERO:
def surreal_ZERO(): return ([],[])

# Definition of surreal ONE:
def surreal_ONE(): return ([surreal_ZERO()],[])

# Definition surreal NEGATIVE ONE:
def surreal_NEGATIVE_ONE(): return ([],[surreal_ZERO()])



######################################################################
# BASIC DEFINITION OF RELATIONSHIPS AND STRUCTURE OF SURREAL NUMBERS:
# 1) surreal_LESS_THAN_OR_EQUAL(X,Y)
# 2) surreal_EQUIVALENT(X,Y)
# 3) surreal_WELL_FORMED_NUMBER(X)
######################################################################

# Definition of the order relationship between two surreal numbers:
def surreal_LESS_THAN_OR_EQUAL(X,Y):
    if X == None or Y == None: return None      # "None" means that there are errors somewhere so STOP.
    (Xs,Xd) = X     # split the first surreal number X into left and right side
    (Ys,Yd) = Y     # split the second surreal number Y into left and right side
    r1 = True
    # if there is a surreal number in the LEFT side of X that is greater than or equal to Y then FALSE
    if Xs != []:
        for Xs_item in Xs:
            if surreal_LESS_THAN_OR_EQUAL(Y,Xs_item):
                r1 = False
                break
    if not r1: return False
    r2 = True
    # if there is a surreal number in the RIGHT side of Y that is less than or equal to X then FALSE
    if Yd != []:
        for Yd_item in Yd:
            if surreal_LESS_THAN_OR_EQUAL(Yd_item,X):
                r2 = False
                break
    return r2
    # Otherwise TRUE

# Definition of the equivalence relationship between two surreal numbers:
def surreal_EQUIVALENT(X,Y):
    # Two surreal numbers X,Y are equivalent if and only if the first is 
    # less than or equal to the second and vice versa:
    return surreal_LESS_THAN_OR_EQUAL(X,Y) and surreal_LESS_THAN_OR_EQUAL(Y,X)

# Check if the surreal number X is well formed:
def surreal_WELL_FORMED_NUMBER(X):
    # A surreal number is well-formed iff:
    # 1) it is a couple
    # 2) the two elements of the pair are two lists
    # 3) every element to the right of X is a well-formed surreal number and similarly it is for the left side
    # 4) no element on the right side of X is less than or equal to an element on the left side.
    if isinstance(X, tuple) and len(X) == 2 and all(isinstance(elem, list) for elem in X):
        (Xs,Xd) = X     # split the surreal number X into left and right side.
        r = True        # recursively evaluate the well-formed property:
        for Xs_item in Xs:
            for Xd_item in Xd:
                if  surreal_WELL_FORMED_NUMBER(Xs_item) and surreal_WELL_FORMED_NUMBER(Xd_item):
                    if surreal_LESS_THAN_OR_EQUAL(Xd_item,Xs_item):
                        print("Bad surreal number format: there are problems, a right member is less than or equal to a left one:")
                        print("LEFT MEMBER -->",Xs_item)
                        print("RIGHT MEMBER -->",Xd_item)
                        quit()
                else:
                        print("Bad surreal number format: one of the two following surreal numbers are not well-formed:")
                        print("LEFT MEMBER -->",Xs_item)
                        print("RIGHT MEMBER -->",Xd_item)
                        quit()
        return True
    else: 
        print("Bad structure of the surreal number:")
        print("-->",X)
        quit()


##################################################################################
# PROCEDURES FOR CONVERTING SURREAL NUMBERS INTO DECIMAL NUMBERS AND VICE VERSA:
# 1)  convert_from_SURREAL_to_decimal_number(X)
# 2)  surreal_DYADIC_FRACTION(n)
# 3)  approximated_inverse_of_decimal_number(x)
# 4)  convert_from_positive_integer_to_surreal_INVERSE(n)
# 5)  surreal_OPPOSITE(X)
# 6)  add_surreal_ONE_to_surreal_X_for_n_times(X,n)
# 7)  surreal_INFINITY()
# 8)  surreal_NEGATIVE_INFINITY()
# 9)  surreal_INFINITESIMAL()
# 10) surreal_NEGATIVE_INFINITESIMAL()
# 11) convert_from_fraction_to_surreal_number(numerator,denominator)
# 12) convert_from_decimal_to_surreal_number(decimal_number)
# 13) reduce_surreal_number(X)
##################################################################################

# Convert from a surreal number X to the equivalent decimal number:
def convert_from_SURREAL_to_decimal_number(X):
    if X == None: return None    # "None" means that there are errors somewhere so STOP.
    if X == surreal_ZERO(): return 0            # Convert surreal ZERO into 0
    if X == surreal_ONE(): return 1             # Convert surreal ONE into 1
    if X == surreal_NEGATIVE_ONE(): return -1   # Convert surreal NEGATIVE ONE into -1
    (LEFTX,RIGHTX) = X          # split the surreal number X into left and right side.
    # Convert into decimal number every surreal number in the left side of X.
    # Get the maximum decimal number from the left side.
    LEFTX_numeric = [convert_from_SURREAL_to_decimal_number(Y) for Y in LEFTX]
    if LEFTX_numeric != []:
        max_LEFTX_number = None
        for leftx_number in LEFTX_numeric:
            if leftx_number == None: return None
            if max_LEFTX_number == None or max_LEFTX_number < leftx_number:
                max_LEFTX_number = leftx_number
    # Convert into decimal number every surreal number in the right side of X.
    # Get the minimum decimal number from the right side.
    RIGHTX_numeric = [convert_from_SURREAL_to_decimal_number(Y) for Y in RIGHTX]
    if RIGHTX_numeric != []:
        min_RIGHTX_number = None
        for rightx_number in RIGHTX_numeric:
            if rightx_number == None: return None
            if min_RIGHTX_number == None or min_RIGHTX_number > rightx_number:
                min_RIGHTX_number = rightx_number
    if LEFTX_numeric == [] and RIGHTX_numeric != []:    # if left side is empty: decrease 1
        return min_RIGHTX_number-1
    elif LEFTX_numeric != [] and RIGHTX_numeric == []:  # if right side is empty: increase 1
        return max_LEFTX_number+1
    elif LEFTX_numeric == [] and RIGHTX_numeric == []:  # if both left and right side are empty: ERROR.
        return None
    else:   # If both left and right sides are not empty: Returns the average of the max and min:
        return (max_LEFTX_number+min_RIGHTX_number)/2   

# It returns the surreal number equivalent to a dyadic fraction 1/(2^n) with n positive integer:
def surreal_DYADIC_FRACTION(n):
    # if n = 0 then 1/2^0 is equal to one so it returns surreal ONE.
    if n == 0: return surreal_ONE()
    # if n > 0 then construct a surreal with the left side ZERO and the right side a dyadic fraction 1/2^(n-1):
    elif n > 0: return ([surreal_ZERO()],[surreal_DYADIC_FRACTION(n-1)])
    # Otherwise: ERROR (if n < 0 then we do not have a fraction but a power of 2).
    else: return None

# If x is a decimal number greather than 1, this procedure returns
# the list of powers of 2 (p1, p2,...pn) for which the sum of the
# related dyadic fractions 1/(2^p1)+1/(2^p2)+...+1/(2^pn) approximates
# the value of 1/x.
# acc returns the approximate value of x obtained.
def approximated_inverse_of_decimal_number(x):
    if x < 1: return None,None
    max_power_of_two = 20   # this is the highest power of 2 considered for the approximation
    inverse_of_x = 1/x      # calculate the inverse of input value x
    list_of_power_of_two = list()   # initializes the list of powers of 2 to be empty
    acc = 0                         # initialise the aproximation value to zero
    for power in range(max_power_of_two):   # for all ther powers:
        dyadic_fraction = 2**(-power)               # calculates every dyadic fraction related to the power
        if acc + dyadic_fraction == inverse_of_x:   # calculate the difference between the aproximation+fraction and x:
            list_of_power_of_two.append(power)      # if aproximation+fraction = x then add to the list the power of 2 and STOP.
            return acc,list_of_power_of_two
        elif acc + dyadic_fraction < inverse_of_x:  # if aproximation+fraction < x then add to the list the power of 2 and CONTINUE.
            list_of_power_of_two.append(power)
            acc += dyadic_fraction
    return acc,list_of_power_of_two                 # if the maximum power of 2 has been reached then STOP.

# It returns true if n is an integer with zero decimal places:
def is_integer(n):
    if isinstance(n, int):
        return True
    if isinstance(n, float):
        return n.is_integer()
    return False

# It returns the surreal number equivalent to the inverse of the positive integer n with a given approximation:
def convert_from_positive_integer_to_surreal_INVERSE(n):
    
    def insert_surreal_DYADIC_FRACTION(L):
        l = len(L)
        if l < 2: return None
        if l == 2:
            return ([surreal_DYADIC_FRACTION(L[1])],[surreal_DYADIC_FRACTION(L[0])])
        elif l > 2:
            return ([insert_surreal_DYADIC_FRACTION(L[1:])],[surreal_DYADIC_FRACTION(L[0])])
        else: return None
        
    # if n is not an integer value or it is less than or equal to 0 then ERROR.
    if not is_integer(n) or n <= 0: return None
    # Get the list of powers of 2 to approximate 1/n with a sum of dyadic fractions:
    _,list_of_power_of_two = approximated_inverse_of_decimal_number(n)
    l = len(list_of_power_of_two)
    # No power of 2 are necessary: it returns surreal ONE:
    if l == 0: return surreal_ONE() 
    # just one power of 2 is enough: it return 1/2^p in surreal form.
    if l == 1: return surreal_DYADIC_FRACTION(list_of_power_of_two[0])
    # If we have a list of many powers of 2: 
    reduced_list_of_powers_of_two = [i for i in list_of_power_of_two]
    reduced_list_of_powers_of_two[l-1] = reduced_list_of_powers_of_two[l-1]-l+1
    for i in range(l-1):
        reduced_list_of_powers_of_two[i] = reduced_list_of_powers_of_two[i]-i-1
    # we can use the reduced list to combine the surreal dyadic fractions:
    return insert_surreal_DYADIC_FRACTION(reduced_list_of_powers_of_two)

# It returns the surreal number equivalent to the opposite of the surreal number X:
def surreal_OPPOSITE(X):
    if X == None: return None           # "None" means that there are errors somewhere so STOP.
    elif X == surreal_ZERO(): return X  # the opposite of surreal ZERO is surreal ZERO.
    elif X == surreal_ONE(): return surreal_NEGATIVE_ONE()  # the opposite of surreal ONE is surreal NEGATIVE ONE
    elif X == surreal_NEGATIVE_ONE(): return surreal_ONE()  # the opposite of surreal NEGATIVE ONE is surreal ONE
    else:
        (LEFTX,RIGHTX) = X  # split the surreal number X into left and right side.
        leftx_result = [surreal_OPPOSITE(x) for x in RIGHTX]    # build the new left side getting the opposite of the X right side.
        rightx_result = [surreal_OPPOSITE(x) for x in LEFTX]    # build the new right side getting the opposite of the X left side.
        return (leftx_result,rightx_result) # combine the new pair and return the new surreal number.

# This function, with X being a surreal number and n being an integer, 
# gets the successor and predecessor operations for surreal numbers.
# If n is a positive integer, it returns the surreal number equivalent to the nth successor of X.
# If n is negative, it returns the surreal number equivalent to the nth predecessor of X.
# If n is 0, it returns X.
def add_surreal_ONE_to_surreal_X_for_n_times(X,n):
    if n == 0: return X # if you add 0 then no change to X
    elif n > 0:         # if n is positive add n to the surreal number X (n-th surreal successor)
        TMP = X
        for i in range(n):
            TMP = ([TMP],[])
        return TMP
    else:               # if n is negative, subtract n from the surreal number X (n-th surreal predecessor)
        TMP = X
        for i in range(-n):
            TMP = ([],[TMP])
        return TMP

# It defines the largest positive number that is considered infinity:
def real_infinity_number(): return 100

# It returns the surreal INFINITY by converting the previous "infinity" number:
def surreal_INFINITY(): return add_surreal_ONE_to_surreal_X_for_n_times(surreal_ZERO(),real_infinity_number())

# It returns the surreal NEGATIVE INFINITY getting the opposite of the previous surreal INFINITY:
def surreal_NEGATIVE_INFINITY(): return surreal_OPPOSITE(surreal_INFINITY())

# It returns the surreal INFINITESIMAL that is a non zero number as smallest as possible based on the biggest number defined:
def surreal_INFINITESIMAL(): return surreal_DYADIC_FRACTION(real_infinity_number())

# It returns the surreal NEGATIVE INFINITESIMAL getting the opposite of the previous surreal INFINITESIMAL:
def surreal_NEGATIVE_INFINITESIMAL(): return surreal_OPPOSITE(surreal_INFINITESIMAL())

# Convert a rational number (both numerator and denominator are integers) into a surreal number with approximation:
def convert_from_fraction_to_surreal_number(numerator,denominator):
    if not is_integer(numerator) or not is_integer(denominator): return None    # Bad input values: ERROR.
    if numerator < 0: 
        if denominator < 0:
            # numerator and denominator are both negative: so the fraction is positive:
            return convert_from_fraction_to_surreal_number(-numerator,-denominator)
        elif denominator == 0:
            # if denominator is 0 and numerator is not 0 then returns surreal NEGATIVE INFINITY.
            return surreal_NEGATIVE_INFINITY()
        else:
            # if numerator is negative and denominator is positive non zero then returns a negative fraction:
            return surreal_OPPOSITE(convert_from_fraction_to_surreal_number(-numerator,denominator))
    elif numerator > 0:
        if denominator < 0:
            # if denominator is negative non zero and numerator is positive then returns a negative fraction:
            return surreal_OPPOSITE(convert_from_fraction_to_surreal_number(numerator,-denominator))
        elif denominator == 0:
            # if denominator is 0 and numerator is not 0 then returns surreal (positive) INFINITY.
            return surreal_INFINITY()
    else:
        if denominator == 0:
            # if both numerator and denominator are zero then ERROR.
            return None
        else:
            # if numerator is zero and denominatore is not zero then returns surreal ZERO.
            return surreal_ZERO()
    # if both numerator and denominator are not zero and positive:
    if numerator > denominator:
        # if the fraction is greater than one extracts the integer value and gets a fraction between 0 and 1:
        new_numerator = numerator % denominator
        number_of_units = numerator // denominator
        # convert to surreal number the fraction between 0 and 1:
        CALCULATION_FRACTION = convert_from_fraction_to_surreal_number(new_numerator,denominator)
        # then add the integer value to the just converted fraction:
        return add_surreal_ONE_to_surreal_X_for_n_times(CALCULATION_FRACTION,number_of_units)
    common_divisor = math.gcd(numerator, denominator)   # here the fraction is necessarely between 0 and 1
    numerator = numerator // common_divisor             # and numerator and denominator are not zero and positive
    denominator = denominator // common_divisor         # we can simplify the fraction.
    # if the fraction is of type 1/m we already have the procedure for the surreal conversion:
    if numerator == 1: return convert_from_positive_integer_to_surreal_INVERSE(denominator)
    # if the fraction is of generic type num/den with num > 1:
    # Algorithm:
    # E.g. 3/5 -> tmp = 3*2 = 6, new_fraction = (6-5)/5 = 1/5 -> ([A_FIFTH],[ONE])
    # E.g. 2/7 -> tmp = 4, new_fraction = -(7-4)/7 = -3/7 -> ([NEGATIVE_THREE_SEVENTHS],[ONE])
    if denominator > 0:
        tmp = numerator * 2
        if tmp > denominator: 
            NEW_FRACTION = convert_from_fraction_to_surreal_number(tmp-denominator,denominator)
            if surreal_LESS_THAN_OR_EQUAL(surreal_ONE(),NEW_FRACTION):
                return ([surreal_ONE()],[NEW_FRACTION])
            else:
                return ([NEW_FRACTION],[surreal_ONE()])
        elif tmp < denominator:
            NEW_FRACTION = surreal_OPPOSITE(convert_from_fraction_to_surreal_number(denominator-tmp,denominator))
            if surreal_LESS_THAN_OR_EQUAL(NEW_FRACTION,surreal_ONE()):
                return ([NEW_FRACTION],[surreal_ONE()])
            else:
                return ([surreal_ONE()],[NEW_FRACTION])
        else:
            return None
    else:
        return None

# Convert a decimal number to a surreal number with approximation:
def convert_from_decimal_to_surreal_number(decimal_number):

    # Convert the decimal number into a fraction with aproximation:
    def convert_from_decimal_to_fraction(decimal_num):
        decimal_precision = 4               # number of reliable digits after the dot
        factor = 10 ** decimal_precision
        rough_numerator = int(round(decimal_num * factor))  # obtain the numerator considering only the reliable digits
        rough_denominator = factor                          # the denominator is just 10^precision 
        common_divisor = math.gcd(rough_numerator, rough_denominator)   # we can simplify the fraction now
        numerator = rough_numerator // common_divisor
        denominator = rough_denominator // common_divisor
        return numerator,denominator
    
    numerator,denominator = convert_from_decimal_to_fraction(decimal_number)    # 1) convert the input decimal number into a fraction
    return convert_from_fraction_to_surreal_number(numerator,denominator)       # 2) convert the obtained fraction into a surreal number

# It reduces the structure and dimension of a surreal number X in a standard and minimum form:
def reduce_surreal_number(X):
    if surreal_WELL_FORMED_NUMBER(X):                       # if X is well-formed
        x = convert_from_SURREAL_to_decimal_number(X)       # we convert it into a decimal number
        return convert_from_decimal_to_surreal_number(x)    # then we re-convert it into a surrel form


####################################################
# FOUR FUNDAMENTAL OPERATIONS WITH SURREAL NUMBERS:
####################################################

# SUM of two surreal numbers X + Y passing through decimal numbers:
def surreal_sum(X,Y):
    if surreal_WELL_FORMED_NUMBER(X) and surreal_WELL_FORMED_NUMBER(Y):
        x = convert_from_SURREAL_to_decimal_number(X)
        y = convert_from_SURREAL_to_decimal_number(Y)
        return convert_from_decimal_to_surreal_number(x+y)

# DIFFERENCE of two surreal numbers X - Y passing through decimal numbers:
def surreal_difference(X,Y):
    if surreal_WELL_FORMED_NUMBER(X) and surreal_WELL_FORMED_NUMBER(Y):
        x = convert_from_SURREAL_to_decimal_number(X)
        y = convert_from_SURREAL_to_decimal_number(Y)
        return convert_from_decimal_to_surreal_number(x-y)    

# PRODUCT of two surreal numbers X * Y passing through decimal numbers:
def surreal_product(X,Y):
    if surreal_WELL_FORMED_NUMBER(X) and surreal_WELL_FORMED_NUMBER(Y):
        x = convert_from_SURREAL_to_decimal_number(X)
        y = convert_from_SURREAL_to_decimal_number(Y)
        return convert_from_decimal_to_surreal_number(x*y)

# DIVISION between two surreal numbers X : Y passing through decimal numbers and infinity:
def surreal_division(X,Y):
    if surreal_WELL_FORMED_NUMBER(X) and surreal_WELL_FORMED_NUMBER(Y):
        x = convert_from_SURREAL_to_decimal_number(X)
        y = convert_from_SURREAL_to_decimal_number(Y)
        X_equivalent_to_ZERO = (x == 0)
        Y_equivalent_to_ZERO = (y == 0)
        if X_equivalent_to_ZERO and Y_equivalent_to_ZERO: return surreal_ONE()
        if X_equivalent_to_ZERO: return surreal_ZERO() 
        if Y_equivalent_to_ZERO:
            if x < 0: return surreal_NEGATIVE_INFINITY()
            else: return surreal_INFINITY()
        return convert_from_decimal_to_surreal_number(x/y)


##################################################
# CALCULATION OF SOME REMARKABLE SURREAL NUMBERS:
##################################################

# Calculate Napier's number in surreal form:
# sum of the series 1/0! + 1/1! + 1/2! + 1/3! + 1/4! + ...
def surreal_E():

    # Calculate the surreal inverse factorial 1/n!
    def inverse_surreal_factorial(n):
        if n == 0: return surreal_ONE()
        elif n > 0: return surreal_product(inverse_surreal_factorial(n-1),convert_from_positive_integer_to_surreal_INVERSE(n))
        else: return None
        
    A = list()
    iterations = 10
    for i in range (iterations):
        A.append(inverse_surreal_factorial(i))
    ACC = surreal_ZERO()
    for X in A:
        ACC = surreal_sum(ACC,X)
    return ACC
    
# Calculate Pi number in surreal form:
# sum of the series 1 - 1/3 + 1/5 -1/7 + 1/9 - 1/11 + ....
# multiplied by 4
def surreal_PI():
    A = list()
    A.append(surreal_ONE())
    for i in range(150):
        if i % 2 != 0:
            A.append(convert_from_positive_integer_to_surreal_INVERSE(3+i*2))
        else:
            A.append(surreal_OPPOSITE(convert_from_positive_integer_to_surreal_INVERSE(3+i*2)))
    ACC = surreal_ZERO()
    for X in A:
        ACC = surreal_sum(ACC,X)
    surreal_FOUR = add_surreal_ONE_to_surreal_X_for_n_times(surreal_ONE(),3)
    return surreal_product(ACC,surreal_FOUR)


####################################################
# CALCULATION OF SOME REMARKABLE SURREAL SEQUENCES:
####################################################

# It calculates a Fibonacci sequence of surreal numbers of lenght n:
# a = ZERO, b = ONE, a = b, b = a+b, recursively
def surreal_fibonacci_sequence(n):
    surreal_sequence = list()
    A, B = surreal_ZERO(), surreal_ONE()
    for _ in range(n):
        surreal_sequence.append(A)
        A, B = B, surreal_sum(A,B)
    return surreal_sequence

# Computes a sequence of surreal pseudo-random numbers of specified length and seed:
# Xn+1 = (Xn * MUL + ADD) always under MAX (by recurring subtractions and not mod())
# where all variables are surreal numbers that are not necessarily integers.
def surreal_random_sequence(seed,length):
    A = convert_from_positive_integer_to_surreal_INVERSE(seed)
    surreal_sequence = list()
    MUL = add_surreal_ONE_to_surreal_X_for_n_times(surreal_ZERO(),3)    # Multiplication factor
    ADD = add_surreal_ONE_to_surreal_X_for_n_times(surreal_ZERO(),7)    # Additive factor
    MAX = add_surreal_ONE_to_surreal_X_for_n_times(surreal_ZERO(),100)  # maximum value
    B = surreal_product(MUL,A)
    C = surreal_sum(B,ADD)
    for i in range(length):
        A = C
        print(i+1,"/",length,")",convert_from_SURREAL_to_decimal_number(A))
        surreal_sequence.append(A)
        B = surreal_product(MUL,A)
        C = surreal_sum(B,ADD)
        while not surreal_LESS_THAN_OR_EQUAL(C,MAX):    # keep the random value under the MAX
            C = surreal_difference(C,MAX)
    return surreal_sequence




###################
###################
## MAIN PROGRAM: ##
###################
###################

if __name__ == "__main__":
    
    # surreal integers from -9 to 9:
    for i in range(-9,10):
        print(i," = ",convert_from_decimal_to_surreal_number(i))
    print()
        
    # Surreal PI:
    PI = surreal_PI()
    print("PI = ",PI,"well formed surreal number:",surreal_WELL_FORMED_NUMBER(PI))
    print("pi = ",convert_from_SURREAL_to_decimal_number(PI))
    print()
    
    # Surreal E:
    E = surreal_E()
    print("E = ",E,"well formed surreal number:",surreal_WELL_FORMED_NUMBER(E))
    print("e = ",convert_from_SURREAL_to_decimal_number(E))
    print()
    
    # 10 elements of the Fibonacci surreal sequence:
    n = 10
    print("Fibonacci sequence (",n,"values):")
    for X in surreal_fibonacci_sequence(n):
        print(convert_from_SURREAL_to_decimal_number(X),X,"well formed surreal number:",surreal_WELL_FORMED_NUMBER(X))
    print()
    
    # Surreal random number sequence:
    n = 10
    seed = 327
    print("Random sequence (",n,"values, seed=",seed,"):")
    for X in surreal_random_sequence(seed,n):
        print(convert_from_SURREAL_to_decimal_number(X),X,"well formed surreal number:",surreal_WELL_FORMED_NUMBER(X))
    
    print("End.")
    






