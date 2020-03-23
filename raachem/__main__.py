from raachem.util.interface import UserInterface
#############################HOOK on ERRORS###############################
import traceback, sys
def show_exception_and_exit(exc_type, exc_value, tb):
	traceback.print_exception(exc_type, exc_value, tb)
	input("Press key to exit.")
	sys.exit(-1)
sys.excepthook = show_exception_and_exit
##########################################################################
UserInterface().render()

