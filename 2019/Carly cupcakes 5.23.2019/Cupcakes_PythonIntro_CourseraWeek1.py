#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:28:04 2019

@author: carly.cheung
"""

#%%
def hello():
    """ prints hello, world """
    print("Hello, world!")

#%%

"""
Exercise:
Write a function 'def areatriangle(b,h)' to compute the area
of a triangle: formula is area = .5*b*h.
Output should look like:
The area of a triangle of base 3 and height 5 is 7.5

You can test your function by executing the following code:
"""
#%%
# The following will test areatriangle()
areatriangle(3,5)
areatriangle(2,20)
#%%
#solution


#%%
To make a string we may use ' or ". Either works equally well. But if the 
string contains one, we need to use the other:
    
#%%

name = "His name is Conan O'Brien"
cat = 'My cat is named "Butters"'
print(name)
print(cat)

#%%
""" 
If you need both a ' and a " in your string, you can use the escape
character \ which tells Python that the following character is to be taken 
as the literal character and is not a quote to delimit the string.  See it
in action escaping the " below:
"""
#%%
both = "My cat's name is \"Butters\""
print(both)

#%%

def fahrenheit_to_celsius(temp):
    """ Converts Fahrenheit temperature to Celsius. 
        Formula is 5/9 of temp minus 32 """
    # Note that this line is not executed
    # end='' keeps print from starting a new line.
    newTemp = 5*(temp-32)/9
    print("The Fahrenheit temperature",temp,"is equivalent to",newTemp,end='')
    print(" degrees Celsius")

#%%
#Test run
fahrenheit_to_celsius(32)
The Fahrenheit temperature 32 is equivalent to 0.0 degrees Celsius

fahrenheit_to_celsius(60)
The Fahrenheit temperature 60 is equivalent to 15.555555555555555 degrees Celsius

    
#%%

"""
Exercise:
Write a function 'def celsius_to_fahrenheit(temp)' to convert Celsius
to Fahrentheit temperature. The formula is (9/5) times temp plus 32. 
Print the output in the form:
The Celsius temperature 50.0 is equivalent to 122.0 degrees Fahrenheit.
"""
#%%
# The following will test the above function
celsius_to_fahrenheit(100)
celsius_to_fahrenheit(0)
celsius_to_fahrenheit(50.)

#%%
#solution

#%%    
def inches_to_feet1(inches):
    """ converts inches to feet and inches """
    feet = inches//12  # division by integer with fraction thrown away
    extra_inches = inches - 12*feet
    print(inches,"inches is",feet,"feet and",extra_inches,"inches") 
#%%
#Test run
inches_to_feet1(52)
52 inches is 4 feet and 4 inches

inches_to_feet1(69)
69 inches is 5 feet and 9 inches

#%%
"""
Exercise: Rewrite inches_to_feet1(inches) calling it inches_to_feet2(inches)
using % to compute the inches. Recall that 19 % 5 will give 4 (the remainder).
Copy and paste the original into the solution area and modify to same typing 
time.
"""
"""
Solution:
"""
#%%

#solution

#%%
    def name():
    """ Input first and last name, combine to one string and print """
    fname = input("Enter your first name: ")
    lname = input("Enter your last name: ")
    fullname = fname + " " + lname

    print("Your name is:", fullname)

#%% 
"""
Exercise:
Extend the name function written in class to include the city and state.
That is, ask two more questions to get the city and the state you live in.
Print where you are from on a new line. Put the customary comma between
city and state. to save time, here is the starting function.
Your run should look like the following (even if this is not the customary
way in your country):
Enter your first name: Bill

Enter your last name: Boyd

Enter the city you live in: Middletown

Enter the state you live in: CT
Your name is: Bill Boyd
You live in:  Middletown, CT

"""
"""
Solution:
"""

#%%
#solution

#%%
def if_statement():
    """ Three slightly difference versions of if: if, if-else, if-elif-else"""
    x = 5
    y = 0
    z = 0
    if x > 0:
        print("x is positive")
        
    if y > 0:
        print("y is positive")
    else:
        print("y is not positive")
        
    # elif can be repeated as often as necessary    
    if z > 0:
        print("z is positive")
    elif z < 0:
        print("z is negative")
    else:
        print("z must be 0")
            
#%%
Exercise:
Write a function absolutevalue(num) that computes the absolute value of
a number. You will need to use an 'if' statement. Remember if a number is 
less than zero then you must multiply by -1 to make it greater than zero.
Give output in the form:

The absolute value of -5  is  5
"""
#%%
# Test runs
absolutevalue(5)
absolutevalue(-5)
absolutevalue(4-4)

#%%
#Solution

#%%
Problem 1_5:
Write a function 'problem1_5(age)'. This function should use if-elif-else
statement to print out "Have a glass of milk." for anyone under 7; "Have
a coke." for anyone under 21, and "Have a martini." for anyone 21 or older.

Tip: Be careful about the ages 7 (a seven year old is not under 7) and 21.
Also be careful to make the phrases exactly as shown for the auto-grader.
"""
#%%
#solution

    
#%%
"""
Test runs (3 of them). Note that the grader program will use different numbers:
problem1_5(5)
Have a glass of milk.

problem1_5(10)
Have a coke.

problem1_5(25)
Have a martini.

#%%
#%%
def fahrenheit_to_celsius1():
    """ BAD. Does not check input before using it. 
    Input from keyboard, which is always a string and must often be
    converted to an int or float. 
    Converts Fahrenheit temp to Celsius.
    """
    
    temp_str = input("Enter a Fahrentheit temperature: ")
    temp = int(temp_str)
    newTemp = 5*(temp-32)/9
    print("The Fahrenheit temperature",temp,"is equivalent to ",end='')
    print(newTemp,"degrees Celsius")

#%%
"""
Test the program above by entering a temperature such as 212. Also check what
happens if you simply press enter.
"""
#%%
def fahrenheit_to_celsius2():
    """ IMPROVED. Does some checking of input before using it.
    Input from keyboard, which is always a string and must often be
    converted to an int or float. 
    Converts Fahrenheit temp to Celsius.
    Uses 'if' to make sure an entry was made.
    """
    
    temp_str = input("Enter a Fahrenheit temperature: ")
    if temp_str:
        temp = int(temp_str)
        newTemp = 5*(temp-32)/9
        print("The Fahrenheit temperature",temp,"is equivalent to ",end='')
        print(newTemp,"degrees Celsius")

#%%
"""
Test the program above by entering the temperature 212 and also by simply
pressing 'Enter' or 'Return' key. Note the improvement. Now try entering 'a'.
"""
#%%
def fahrenheit_to_celsius3():
    """ MORE IMPROVED. Does even more checking of input before using it. 
    Input from keyboard, which is always a string and must often be
    converted to an int or float. 
    Converts Fahrenheit temp to Celsius.
    Uses if to check whether input is a number and then uses .isdigit() method 
    of strings to check whether input is made of of digits. 
    """
        
    temp_str = input("Enter a Fahrentheit temperature: ")
    if temp_str:
        if temp_str.isdigit():  
            temp = int(temp_str)
            newTemp = 5*(temp-32)/9
            print("The Fahrenheit temperature",temp,"is equivalent to ",end='')
            print(newTemp,"degrees Celsius")
        else:
            print("You must enter a number. Bye")
#%%
"""
Test the program above by entering the temperature 212, by simply pressing 
'Enter' or 'Return' key, and by entering 'a'. Note the improvement. We will
leave the function at this point though further improvements could be made.
"""
#%%
Problem 1_7:
Write a function problem1_7() that computes the area of a trapezoid. Here is the
formula: A = (1/2)(b1+b2)h. In the formula b1 is the length of one of the 
bases, b2 the other. The height is h and the area is A. Basically, this 
takes the average of the two bases times the height. For a rectangle b1 = b2, 
so this reduces to b1*h. This means that you can do a pretty good test of the 
correctness of your function using a rectangle (that way you can compute the 
answer in your head). Use input statements to ask for the bases and the height.
Convert these input strings to real numbers using float(). Print the output
nicely EXACTLY like mine below.

Tip: Be careful that your output on the test case below is exactly as shown
so that the auto-grader judges your output correctly.  The auto-grader does
not look at your input statements, so you don't have to use my input prompts
if you don't want to. However, the auto-grader will enter the three inputs in
the order shown. See the other test run below.

problem1_7()

Enter the length of one of the bases: 3

Enter the length of the other base: 4

Enter the height: 8
The area of a trapezoid with bases 3.0 and 4.0 and height 8.0 is 28.0

"""  
#%%
def problem1_7():
    b11 = input("Enter the length of one of the bases:")
    if b11:
        if b11.isdigit():
            b1 = int(b11)
        else: 
            print("You must enter a number.")
    b22 = input("Enter the length of the other base:")
    if b22:
        if b22.isdigit():
            b2 = int(b22)
        else: 
            print("You must enter a number.")
    hh = input("Enter the height:")
    if hh:
        if hh.isdigit():
            h = int(hh)
    else: 
        print("You must enter a number.")
    A = (1/2)*(b1+b2)*h
    print("The area of a trapezoid with bases", float(b1), "and", float(b2), "and height", float(h), "is", A)

#%%
def cheer():
    """ Prints 2 4 6 8, who do we appreciate .... Note that everything in 
    the while loop is indented. The first line not indented is the first
    line following the while loop. """
    ct = 2
    while ct <= 8:
        print(ct,end=" ")  # end = " " keeps from starting a new line
        ct = ct + 2
    print()                # now we'll start a new line
    print("Who do we appreciate?")
    print("COURSERA!")
 
#%%   
"""
Exercise:
Write a function count_down() that starts at 10 and counts down to rocket 
launch. It's output should be 10 9 8 7 6 5 4 3 2 1 BLASTOFF! You can make
all the numbers on the same line or different lines. Use a while loop.
"""
#%%
#solution

#%%
The 'for' loop. This loop uses an iterator to determine how many times to go
through the loop. The iterator we use below is 'range(start, stop, step)'.
"""
#%%
def cheer2():
    """ Same as cheer, but uses a for loop and range()
    range uses a start number, a stop number and a step size. """
    for ct in range(2,9,2):
        print(ct,end=' ')
    print()
    print("Who do we appreciate?")
    print("COURSERA!")
  
#%%  
"""
Exercise:
Write a function countdown1() that starts at 10 and counts down to rocket 
launch. It's output should be 10 9 8 7 6 5 4 3 2 1 BLASTOFF! You can make
all the numbers on the same line or different lines. Use a 'for' loop and
range(). range has a start and a stop and a step that MAY BE NEGATIVE.

"""
#%%
#solution

#%%
Problem 1_3:  
Write a function problem1_3(n) that adds up the numbers 1 through n and
prints out the result. You should use either a 'while' loop or a 'for' loop.
Be sure that you check your answer on several numbers n.  Be careful that your
loop steps through all the numbers from 1 through and including n.

Tip: As this involves a loop you could make an error that causes it to run 
forever. Usually Control-C will stop it. See suggestions at the beginning of 
this document.  With loops take care that your first and last iterations are
what you expect. A print statement can be inserted in the loop to monitor it, 
but be sure this isn't in the submitted function.
"""
#%%


    

#%%
def problem1_3(n):
    my_sum = 0
    for i in range (n+1):
        my_sum = i + my_sum 
    print(my_sum)

#%%

"""
Test run. Note that the grader program will use a different number for n:

problem1_3(6)
21

#%%