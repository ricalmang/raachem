#!/usr/bin/env python
import os
current_folder=str(os.path.dirname(os.path.abspath(__file__)))

def evaluate_list(folder,opt_a,scripts=[]):
	exetensions = [".xyz",".chk",".GAUSSIAN"]
	a_options = ["2","3","4"]
	for file in os.listdir(folder):
		if os.path.isdir(os.path.join(folder,file)):
			if os.path.join(folder,file).endswith("miniconda3"):continue
			if os.path.join(folder,file).endswith("bin"):continue
			evaluate_list(os.path.join(folder,file),opt_a,scripts)
			continue
		try:
			ext = os.path.splitext(file)[1]
			if any(ext == a for a,b in zip(exetensions,a_options) if b in opt_a):
				print(os.path.join(folder,file))
				scripts.append(os.path.join(folder,file))
			if "1" not in opt_a: continue
			if len(ext) != 7: continue
			if not ext[1].isalpha():continue
			if not ext[1].islower():continue
			if not all(a.isdigit() for a in ext[2:]):continue
			else:print(os.path.join(folder,file)) 
			scripts.append(os.path.join(folder,file))
		except IndexError:
			print(file," has no extension!")
			continue
	return	scripts
print("What file types do you want to remove based on their extension?")
print("0 - To return without delleting files")
print("1 - Files with the '.XYYYYY' extension where X are letters and Y are numbers")
print("2 - Files with the '.xyz' extension")
print("3 - Files with the '.chk' extension")
print("4 - Files with the '.GAUSSIAN' extension")
while True:
	opt_a = input().split()
	if all(a in ["0","1","2","3","4"] for a in opt_a):break
	else: print("Invalid Input, please try again")
if "0" in opt_a: exit()
print("Files to be removed:")
filename = evaluate_list(current_folder, opt_a)
print("Do you want to remove the above files?(y/n)")
option = input()
if option == "y":
	for file in filename:
		os.system("rm {0}".format(file))
	print("Done!")
else:
	print("No files will be deleted")

