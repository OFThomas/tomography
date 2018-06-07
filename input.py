#############################################################
# Title: User input
#
# Date created: 7th June 2018
#
# Language:    Python 3
#
# Overview:    
#
# Details:
#
# Usage: to be imported
#
#############################################################

# Ask the user a yes/no question
def yes_or_no(string):
    response = input("\n" + string + " [y/n]")
    while 1==1:
        if response == 'y' or response == "Y":
            print("YES")
            return 0
        elif response == 'n' or response == "N":
            print("NO")
            return 1
        else:
            print(response,"isn't a valid response; try again... ")
            response = input("\n" + string + " [y/n]")

# Get an integer from the user
#   NEEDS TESTING!!
#   Argument 1: A message to print for the benefit of the user
#   Argument 2: The type of the requested value
def get_user_value(string, type):
    while 1==1:
        response = input(string + " " + "[" + type + "]:")
        try:
            if type == "integer" : user_val = int(response)
            elif type == "float" : user_val = float(response)
            elif type == "complex" : user_val = complex(response)
            return user_val
        except ValueError:
            print(response,"is not of type","[" + type + "]:","; try again...")

# Get float from the user
# NEEDS TESTING!!
def get_user_float(string):
    return float(input(string + " "))