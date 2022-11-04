from functions import *

# The algorithm is divided into two steps depending on the parity of the order of the tree
n = ["even","odd"] 

begin = int(input("Initial value for k: "))
end = int(input("Final value for k: "))

file = open(f"k: {begin} -- {end}.txt", "w")
file.write(f"---------------- RESULTS OF THE TEST TO K IN [{begin}, {end}] ----------------\n\n\n")

for k in range(begin, end + 1):
    file.write(f"{dots}k = {k}:\n{dots}\n")

    for elem in n:
        # The function sets the position i of the last vertex q+i to be analyzed
        limit = setLimit(elem, k)
        # The function sets the position i of the first vertex q+i to be analyzed
        initialValue = setInitialValue(elem)
        
        for i in range(initialValue, limit + 1):
            if(i == initialValue):
                if(elem == "even"):
                    file.write(f"Order: {elem} (2q)\n")
                else:
                    file.write(f"Order: {elem} (2q + 1)\n")
            # The function calculates the polynomial whose roots belonging to the interval (0,1) are 1-a(T)
            p = calcPol(elem, i, k)
            
            # The function selects the valid roots of the polynomial p and converts them to the respective a(T)
            list_aT = calcAlgebraicConnectivity(p)
            
            file.write(f"\tVertex q + {i}:\n")
            file.write(f"\t\tPolynomial: {p}\n")
            file.write(f"\t\tCalculated a(T): {list_aT}\n")

            # The function checks if there is any calculated algebraic connectivity 
            # whose difference by theoretical algebraic connectivity is less than the tolerance
            verifyAlgebraicConnectivity(file, (n, k, i), list_aT)

file.close()