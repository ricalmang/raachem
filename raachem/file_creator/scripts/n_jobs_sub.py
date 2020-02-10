#!/usr/bin/env python
import os
directory=os.path.dirname(os.path.abspath(__file__))

def file_weeder(ext_to_weed):
	"""Looks up files with the extensions provided in current directory"""
	fulllist=os.listdir(os.path.dirname(os.path.abspath(__file__)))
	weeded_list=[]
	for extension in ext_to_weed:
		matching = [s for s in fulllist if extension in s]
		for match in matching:
			weeded_list.append(match)
	return sorted(weeded_list)

	
anr_script , files_for_p = [], []
all_files = file_weeder([".pbs"])
print("The following files will compose the script")
for idx,i in enumerate(all_files):
	if idx == 0:
		counter , a = 0, [] 
	counter += 1
	a.append(i)
	if counter == 3 or idx + 1 == len(all_files):
		files_for_p.append(a)
		counter , a = 0, [] 	
for line in files_for_p:
	print(" ".join(["{0}".format(a) for a in line]))

while True:
	n_items = input("How many of these items do you want to submit?\n 'a' for all\n'0' to cancel\n")
	if n_items == "a":
		n_items = int(len(all_files))
		break
	elif n_items == "0":
		print("Ok!, Exiting now!")
		exit()
	elif any([n_items == str(i+1) for i in range(len(all_files))]):
		n_items = int(n_items)
		break

for file in file_weeder([".pbs"])[:n_items]:
	anr_script.append("qsub {0}".format(file))
	anr_script.append("sleep 0.2")
	anr_script.append("rm {0}".format(file))
	anr_script.append("sleep 0.2")
anr_script.append("rm sub_s_name")

if anr_script:
	with open(os.path.join(directory,"sub_s_name"),"w+") as script:
		script.write("\n".join(anr_script))
os.system("chmod 755 sub_s_name")


