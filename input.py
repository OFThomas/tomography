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
def get_user_int(string):
    return int(input(string + " "))
