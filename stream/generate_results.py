import subprocess

cmd_list = ["make all"]

file_name = ["serial_stream", "parallel_stream", "parallel_db_edt_breakdown_stream"]
#file_name = ["native_stream", "serial_stream", "sertial_edt_breakdown_stream", "parallel_stream", "parallel_db_breakdown_stream", "parallel_db_edt_breakdown_stream"]

db_size = [i*10**j for i in range(1, 10) for j in range(5, 7)]
# <file_name> <db_size> <split_size>
#export = [[j + str(i) for i in db_size] for j in file_name]
#export = ["ns", "ss", "sebs", "ps", "pdbbs", "pdbebs"]
export_path = "./results/"
export = [export_path + i for i in ["ss.dat", "ps.dat", "pdbebs.dat"]]
iteration = [5]
split_size = [10**0, 10**1, 10**2, 10**3]
scalar = 3.0
verify = True
verbose = False
 
trial = 5

#ss = ["./serial_stream -d 10000-i -e ss"]
# "./" + <file_name> + -d db_size + -e <export_name> -i <iteration> -p <split_size> -verify <flag> -s scalar <value> -v <flag>
#ss = ["./" + f + " -d " + str(d) + " -e " + e + " -i " + str(i) + " -p " + str(p) + " -s " + str(s)  for f in file_name for d in db_size for e in export for i in iteration for p in split_size for s in scalar]

ss = ["./" + file_name[0] + " -d " + str(d) + " -e " + export[0] + " -i " + str(i) for d in db_size for i in iteration]
ss = [i + " -s " + str(scalar) for i in ss] if scalar != 3.0 else ss
ss = [i + " -r" for i in ss] if verify == True else ss
ss = [i + " -v" for i in ss] if verbose == True else ss

#print ss

ps = ["./" + file_name[1] + " -d " + str(d) + " -e " + export[1] + " -i " + str(i) for d in db_size for i in iteration]
ps = [i + " -s " + str(scalar) for i in ps] if scalar != 3.0 else ps
ps = [i + " -r" for i in ps] if verify == True else ps
ps = [i + " -v" for i in ps] if verbose == True else ps

#print ps

pdbebs = ["./" + file_name[2] + " -d " + str(d) + " -e " + export[2] + " -i " + str(i) for d in db_size for i in iteration]
pdbebs = [i + " -p " + str(10) for i in pdbebs]
pdbebs = [i + " -s " + str(scalar) for i in pdbebs] if scalar != 3.0 else pdbebs
pdbebs = [i + " -r" for i in pdbebs] if verify == True else pdbebs
pdbebs = [i + " -v" for i in pdbebs] if verbose == True else pdbebs


print pdbebs
# run each implementation (running  iteration times) trial times
[subprocess.call(i.split(" ")) for i in pdbebs]

plot_cmd = [(file_name[i], "gnuplot  -e  set terminal pngcairo; set output \'" + file_name[i] + ".png\'; plot \'" + export[i] + "\' using 1:2") for i in range(len(file_name))]


selected = ["parallel_db_edt_breakdown_stream"]
[subprocess.call(j.split("  ")) for i, j in plot_cmd if i in selected]
#print [j.split("  ") for i, j in plot_cmd if i in selected]


# Need to take care of split_size

#Need to plot Native vs serial vs parallel vs parallel_db_edt_breakdown
