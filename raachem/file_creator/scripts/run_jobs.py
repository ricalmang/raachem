#!/usr/bin/env python
import os
current_folder=str(os.path.dirname(os.path.abspath(__file__)))
print("========================")
print("====Enter 0 to cancel===")
print("========================")
def evaluate_list(folder,level=0,scripts=[" "]):
	for file in os.listdir(folder):
		if os.path.isdir(os.path.join(folder,file)):
			#print("{0:}{1:}".format(" "*level, file))
			level += 1
			evaluate_list(os.path.join(folder,file),level,scripts)
		if file == "sub_s_name":
			scripts.append(os.path.join(folder,file))
			print("{0:<2}{1:}".format(len(scripts)-1," - " + scripts[-1]))
			os.system("chmod 755 {0:}".format(os.path.join(folder,file)))
	return	scripts
try:	
	filename = evaluate_list(current_folder)[int(input())]
	if not filename == " ":
		os.chdir(filename.replace("sub_s_name",""))
		os.system("{0:}".format(filename))
		#os.rename(filename,filename.replace("_sub","_subed"))
	print("Done!")

except IndexError:
	print("Not a valid entry!")
except NameError:
	print("Should be a number!")
