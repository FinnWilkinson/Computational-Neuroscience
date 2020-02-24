#for the submission uncomment the submission statements
#so submission.README

from math import *
import numpy as np

from submission import *

def seven_segment(pattern):

    def to_bool(a):
        if a==1:
            return True
        return False


    def hor(d):
        if d:
            print(" _ ")
        else:
            print("   ")

    def vert(d1,d2,d3):
        word=""

        if d1:
            word="|"
        else:
            word=" "

        if d3:
            word+="_"
        else:
            word+=" "

        if d2:
            word+="|"
        else:
            word+=" "

        print(word)



    pattern_b=list(map(to_bool,pattern))

    hor(pattern_b[0])
    vert(pattern_b[1],pattern_b[2],pattern_b[3])
    vert(pattern_b[4],pattern_b[5],pattern_b[6])

    number=0
    for i in range(0,4):
        if pattern_b[7+i]:
            number+=pow(2,i)
    print(int(number))

submission=Submission("Finn_Wilkinson")
submission.header("Finn Wilkinson")

six=[1,1,-1,1,1,1,1,-1,1,1,-1]
three=[1,-1,1,1,-1,1,1,1,1,-1,-1]
one=[-1,-1,1,-1,-1,1,-1,1,-1,-1,-1]

seven_segment(three)
seven_segment(six)
seven_segment(one)

weight_matrix = np.asarray([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                            [0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                            [0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                            [0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                            [0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                            [0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                            [0.0,0.0,0.0,0.0,0.0,0.0,0.0]])

for i in range(0,7) :
    for j in range(0,7) :
        if (i != j) :
            weight_matrix[i,j] = (1/3) * ((six[i]*six[j]) + (three[i]*three[j]) + one[i]*one[j])

print(weight_matrix)

##this assumes you have called your weight matrix "weight_matrix"
submission.section("Weight matrix")
submission.matrix_print("W",weight_matrix)

def McCulloch_Pitts_Formula(weight_matrix, test) :
    while (1) :
        next = [0.0,0.0,0.0,0.0,0.0,0.0,0.0, test[7], test[8], test[9], test[10]]
        for i in range(0,7) :
            for j in range(0,7) :
                next[i] = next[i] + (weight_matrix[i,j] * test[j] )

            if (next[i] > 0) :
                next[i] = 1

            else :
                next[i] = -1

        if (next == test) :
            seven_segment(test)
            submission.seven_segment(test)
            return next

        else :
            test = next
        seven_segment(test)
        submission.seven_segment(test)


print("test1")
submission.section("Test 1")

test=[1,-1,1,1,-1,1,1,-1,-1,-1,-1]

seven_segment(test)
submission.seven_segment(test)
test = McCulloch_Pitts_Formula(weight_matrix, test)

##for COMSM0027

##where energy is the energy of test
#submission.print_number(energy)

##this prints a space
#submission.qquad()

#here the network should run printing at each step
#for the final submission it should also output to submission on each step

print("test2")

test=[1,1,1,1,1,1,1,-1,-1,-1,-1]
submission.section("Test 1")

seven_segment(test)
submission.seven_segment(test)
test = McCulloch_Pitts_Formula(weight_matrix, test)

##for COMSM0027
##where energy is the energy of test
#submission.print_number(energy)

##this prints a space
#submission.qquad()

#here the network should run printing at each step
#for the final submission it should also output to submission on each step


submission.bottomer()
