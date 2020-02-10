import os, math
from difflib import SequenceMatcher
from raachem.util.constants import keywords, cf
from raachem.file_class.gjf import GjfFile
from raachem.util.gen_purp import file_weeder, w_any, read_item, Var, sel_files

def pbs_creator():
	def validate_gjf(weeded_list):
		heavy_e = Var().heavy_atom
		if True:
			print("---------------------------------------------------------------------------")
			print("{:^30}{:^15}{:^10}{:^10}{:^10}".format("File","e- number","charge","multip","Validated"))
			print("---------------------------------------------------------------------------")
		novel_keys = []
		for item in weeded_list:
			gjf = GjfFile(read_item(item))
			split_list = [i for i in gjf.list[1:gjf.title_idx()] if not i.lower().startswith("%chk")]
			for x in [None,"/","(",")",",","=","%",":"]:
				split_list = [a for b in [i.split(x) for i in split_list] for a in b if len(a) > 3]
			no_match = [i for i in split_list if i.lower() not in [j.lower() for j in keywords] and not i[0].isdigit()]
			typo_error = []
			for b in no_match:
				p_matches = [[b,i,SequenceMatcher(None,b,i).ratio()] for i in keywords if SequenceMatcher(None,b,i).ratio() > 0.6]
				if p_matches:
					p_matches.sort(key=lambda x: x[2])
					typo_error.append([p_matches[0][0],p_matches[0][1]])
			novel_keys.append([i for i in no_match if i not in [a[0] for a in typo_error]])
			print(" {:<30}{:>10}{:>10}{:>10}{:^15}".format(gjf.name(),gjf.n_electrons(),gjf.charge(),gjf.multiplicity(),gjf.c_m_validate_txt()))
			if any([gjf.basis_errors(),typo_error,gjf.ecp_errors(heavy_e),len(gjf.name().split()) != 1]):
				print("{:>7}+----------------------------ALERT---------------------------+".format(" "))
				for error in gjf.basis_errors(): print("{:>8}{:>60}".format("|",error+" |"))
				for error in gjf.ecp_errors(heavy_e): print("{:>8}{:>60}".format("|",error+" |"))
				for i in typo_error: print("{:>8}{:>60}".format("|", "Is '{}' a typo of '{}'?".format(i[0],i[1]) + " |"))
				if len(gjf.name().split()) != 1: print("{:>8}{:>60}".format("|", "Filename must not contain spaces!" + " |"))
				print("{:>7}+-----------------------------------------------------------+".format(" "))
		print("---------------------------------------------------------------------------\n")
		novel_keys = list(dict.fromkeys([a for b in novel_keys for a in b]))
		if novel_keys:
			print("The following keys were not recognized:")
			print(novel_keys)
			print("---------------------------------------------------------------------------\n")
	def pbs_creator_athene(weeded_list):
		outscript=[]
		outscript.append("qstat \nsleep 1 \n")
		for i in weeded_list:
			gjf = GjfFile(read_item(i))
			if  gjf.n_proc() == None:
				print("ERROR ON FILE:", i, ", number of processors could not be identified!")
				continue
			if  gjf.mem() == None:
				print("ERROR ON FILE:", i, ", memory request could not be identified!")
				continue
			else:
				out2bi=[]
				out2bi.append("#!/bin/bash\n")
				out2bi.append("##")
				out2bi.append("## queue/job setup, adapt these to your needs!")
				out2bi.append("##\n")
				out2bi.append("                                        # desired queue")
				out2bi.append("#PBS -q rehbein")
				out2bi.append("                                        # number of nodes, cores per node")
				out2bi.append("#PBS -l nodes=1:ppn="+str(gjf.n_proc()))
				out2bi.append("                                        # estimated runtime")
				out2bi.append("#PBS -l walltime=24:00:00")
				out2bi.append("                                        # max. memory per process")
				out2bi.append("#PBS -l pmem="+str(math.ceil(gjf.mem()/gjf.n_proc()))+"mb\n")
				out2bi.append("#ATHOS -o req-cpuset\n")
				out2bi.append("#ATHOS -deps athos-gaussian\n")
				out2bi.append("source /software/g09.ifort.avx.mkl/bsd/g09.profile")
				out2bi.append("PATH=$PATH:/software/g09.ifort.avx.mkl:/software/NBO6/bin")
				out2bi.append("export PATH\n")
				out2bi.append("cd /home/$USER{}".format(str(os.path.basename(cf))))
				out2bi.append("g09 "+ i + " > " + i.replace(".gjf",".log") + "           # start Gaussian09\n")
				w_any(out2bi,"w",i.replace(".gjf",".pbs"))
				outscript.append("qsub {0} \nsleep 1 \nrm {0} \nsleep 1 ".format(i.replace(".gjf",".pbs")))
				print("PBS script created for: {}".format(i))
		outscript.append("qstat")
		w_any(outscript,"w",variables.sub_s_name)
		print("\nJob done!\n")
		return
	def pbs_creator_usp(weeded_list):
		outscript = []
		outscript.append("qstat -u $USER\nsleep 1 \n")
		while True:
			print("Which queue?")
			print("0 - Cancel")
			print("1 - Short")
			print("2 - Default")
			print("3 - Held Default")
			print("4 - Sequential Default")
			queue = [input()]
			if queue[0] in [str(a) for a in range(5)]: break
			else: print("\nCould not understand request: "+ queue[0])
		if queue[0] == "0":	return
		elif queue[0] == "1": queue.extend(["short", "24", ""])
		elif queue[0] == "2": queue.extend(["default", "300", ""])
		elif queue[0] == "3": queue.extend(["default", "300", "-h"])
		elif queue[0] == "4": queue.extend(["default", "300", None])

		for idx,i in enumerate(weeded_list):
			gjf = GjfFile(read_item(i))
			if gjf.n_proc() == None:
				print("ERROR ON FILE:"+i+ ", number of processors could not be identified!")
				continue
			else:
				out2bi=[]
				out2bi.append("#!/bin/bash")
				out2bi.append("#PBS -N "+ i.replace(".gjf",""))
				out2bi.append("#PBS -l nodes=1:ppn=" + str(gjf.n_proc()))
				out2bi.append("#PBS -q {}".format(queue[1]))
				out2bi.append("#PBS -l walltime={}:00:00".format(queue[2]))
				out2bi.append("#PBS -m ae")
				if variables.heimdall_notification: out2bi.append("#PBS -M {}".format(variables.heimdall_mail))
				out2bi.append('echo -e "\\n## Job iniciado em $(date +%d-%m-%Y_%T) #####################\\n"')
				out2bi.append('INP=$PBS_JOBNAME".gjf"')
				out2bi.append('OUT=$PBS_JOBNAME".log"')
				out2bi.append('echo -e "\\n## Jobs ativos de $USER: \\n"')
				out2bi.append("qstat -an -u $USER")
				out2bi.append('echo -e "\\n## Node de execucao do job:         $(hostname -s) \\n"')
				out2bi.append('echo -e "\\n## Numero de tarefas para este job: $PBS_TASKNUM \\n"')
				out2bi.append("module load softwares/gaussian/09d01-pgi-14.3")
				out2bi.append("cd $PBS_O_WORKDIR")
				out2bi.append("EXE=g09")
				out2bi.append('echo -e "\\n## Diretorio de submissao do job:   $PBS_O_WORKDIR \\n"')
				out2bi.append('echo -e "\\n## Diretorio de scratch do job:     $GAUSS_SCRDIR \\n"')
				out2bi.append('echo -e "\\n## Arquivo de input:                $INP \\n"')
				out2bi.append("time $EXE < $INP > $OUT")
				out2bi.append('echo -e "\\n## Job finalizado em $(date +%d-%m-%Y_%T) ###################"')
				out2bi.append('qstat -f $PBS_JOBID')
				out2bi.append(' ')
				w_any(out2bi,"w",i.replace(".gjf",".pbs"))
				if queue[0] in ["1","2","3"]:
					outscript.append("qsub {1} {0}  \nsleep 0.2 \nrm {1} \nsleep 0.2\n".format(queue[-1], i.replace(".gjf",".pbs")))
				elif queue[0] in "4":
					if idx == 0:
						outscript.append("job_{0}=$(qsub {1}) \necho $job_{0}\nsleep 0.2 \n".format(idx,i.replace(".gjf",".pbs")))
					elif idx > 0:
						outscript.append("job_{0}=$(qsub -W depend=afterany:$job_{2} {1}) \necho $job_{0}\nsleep 0.2 \n".format(idx,i.replace(".gjf",".pbs"),idx-1))
				print("PBS script created for: {}".format(i))
		outscript.append("qstat -u $USER\n")
		outscript.append("rm {}".format(variables.sub_s_name))
		w_any(outscript,"w","{}".format(variables.sub_s_name))
		print("\nJob done!\n")
	def pbs_creator_lcca(weeded_list):
		outscript =	[]
		outscript.append("squeue -u $USER\nsleep 1 \n")
		while True:
			print("Which queue?")
			print("0 - Cancel")
			print("1 - Short")
			print("2 - Default")
			queue=[input()]
			if queue[0] in [str(a) for a in range(3)]: break
			else: print("\nCould not understand request: " + queue[0])
		if queue[0] == "0": return
		elif queue[0] == "1": queue.extend(["short","24"])
		elif queue[0] == "2": queue.extend(["default","100"])

		for i in weeded_list:
			gjf = GjfFile(read_item(i))
			if  gjf.n_proc() == None:
				print("ERROR ON FILE:"+i+ ", number of processors could not be identified!")
				continue
			out2bi=[]
			out2bi.append("#!/bin/bash -f")
			out2bi.append("#SBATCH --partition=SP2")
			out2bi.append("#SBATCH --ntasks=1")
			out2bi.append("#SBATCH --cpus-per-task={}".format(gjf.n_proc()))
			out2bi.append("#SBATCH -J {}".format(i.replace(".gjf","")))
			out2bi.append("#SBATCH --time={}:00:00".format(queue[2]))
			out2bi.append(" ")
			out2bi.append("#OpenMP settings:")
			out2bi.append("export OMP_NUM_THREADS=1")
			out2bi.append("export MKL_NUM_THREADS=1")
			out2bi.append("export OMP_PLACES=threads")
			out2bi.append("export OMP_PROC_BIND=spread")
			out2bi.append(" ")
			out2bi.append("echo $SLURM_JOB_ID      ")
			out2bi.append("echo $SLURM_SUBMIT_DIR")
			out2bi.append("echo $SLURM_JOB_NODELIST")
			out2bi.append("echo $SLURM_NTASKS")
			out2bi.append(" ")
			out2bi.append("export PATH=/scratch/apps/pgi/linux86-64/18.10/bin:$PATH")
			out2bi.append("export LD_LIBRARY_PATH=/scratch/apps/pgi/linux86-64/18.10/lib:$LD_LIBRARY_PATH")
			out2bi.append("export INCLUDE=/scratch/apps/pgi/linux86-64/18.10/include:$INCLUDE")
			out2bi.append("export PATH=/scratch/apps/pgi/g09:$PATH")
			out2bi.append("export LD_LIBRARY_PATH=/scratch/apps/pgi/g09:$LD_LIBRARY_PATH")
			out2bi.append("mkdir -p /scratch/$USER/g09/workdir")
			out2bi.append("export GAUSS_SCRDIR=/scratch/$USER/g09/workdir")
			out2bi.append("export OMP_PROC_BIND=FALSE")
			out2bi.append("export g09root=/scratch/apps/pgi")
			out2bi.append("source /scratch/apps/pgi/g09/bsd/g09.profile")
			out2bi.append("")
			out2bi.append("cd $SLURM_SUBMIT_DIR")
			out2bi.append("g09 < ./{0}.gjf > {0}.log".format(i.replace(".gjf","")))
			out2bi.append(' ')
			w_any(out2bi,"w",i.replace(".gjf",".sh"))
			outscript.append("sbatch {0} \nsleep 0.2 \nrm {0}\nsleep 0.2\n".format(i.replace(".gjf",".sh")))
			print("PBS script created for: {}".format(i))
		outscript.append("squeue -u $USER")
		w_any(outscript,"w","{}".format(variables.sub_s_name))
		print("\nJob done!\n")
	def option():
		input("Enter to confirm")
		files = file_weeder([".gjf"]) if variables.folder_op else sel_files(file_weeder([".gjf"]))
		if not files: return
		print("0 - Go back")
		print("1 - For Athene server")
		print("2 - For USP server")
		print("3 - For LCCA server")
		option=input()
		if option == "0": return
		elif option == "1": pbs_creator_athene(files)
		elif option == "2": pbs_creator_usp(files)
		elif option == "3": pbs_creator_lcca(files)
		else: print("Invalid input. Could not process request!")
	validate_gjf(file_weeder([".gjf"]))
	variables = Var()
	option()


