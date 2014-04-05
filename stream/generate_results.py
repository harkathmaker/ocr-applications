import subprocess

subprocess.call("make all".split(" "))

file_name = ["serial_stream", "parallel_stream", "parallel_db_edt_breakdown_stream"]
#file_name = ["native_stream", "serial_stream", "sertial_edt_breakdown_stream", "parallel_stream", "parallel_db_breakdown_stream", "parallel_db_edt_breakdown_stream"]

db_size = [i*10**j for i in range(1, 10) for j in range(5, 7)]
#db_size = [100]
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

# "./" + <file_name> + -d db_size + -e <export_name> -i <iteration> -p <split_size> -r -s scalar <value> -v
ss = ["./" + file_name[0] + " -d " + str(d) + " -e " + export[0] + " -i " + str(i) for d in db_size for i in iteration]
ss = [i + " -s " + str(scalar) for i in ss] if scalar != 3.0 else ss
ss = [i + " -r" for i in ss] if verify == True else ss
ss = [i + " -v" for i in ss] if verbose == True else ss
#print ss

ps = ["./" + file_name[1] + " -d " + str(d) + " -e " + export[1] + " -i " + str(i) for d in db_size for i in iteration]
ps = [i + " -p " + str(10) for i in ps]
ps = [i + " -s " + str(scalar) for i in ps] if scalar != 3.0 else ps
ps = [i + " -r" for i in ps] if verify == True else ps
ps = [i + " -v" for i in ps] if verbose == True else ps
#print ps

pdbebs = ["./" + file_name[2] + " -d " + str(d) + " -e " + export[2] + " -i " + str(i) for d in db_size for i in iteration]
pdbebs = [i + " -p " + str(10) for i in pdbebs]
pdbebs = [i + " -s " + str(scalar) for i in pdbebs] if scalar != 3.0 else pdbebs
pdbebs = [i + " -r" for i in pdbebs] if verify == True else pdbebs
pdbebs = [i + " -v" for i in pdbebs] if verbose == True else pdbebs
#print pdbebs

# run each implementation (running  iteration times) trial times
[subprocess.call(i.split(" ")) for i in ss + ps + pdbebs]

#plot_prefix = "gnuplot  -e  set terminal pngcairo; set output \'"
#plot_dest = "./plots/"
#title = "set title \"Datablock Sizes vs Performance Time\";"
#x_label = "set xlabel \"Datablock Size\";"
#y_label = "set ylabel \"Time (sec)\";"
#x_range = "set xrange [0:100];" 
#y_range = "set y range
#plot_cmd = [(file_name[i], plot_prefix + plot_dest + file_name[i] + ".png\';" + title + x_label + y_label + " plot \'" + export[i] + "\' using 1:2") for i in range(len(file_name))]
#plot_cmd = [(file_name[i], "gnuplot  -e  set terminal pngcairo; set output \'" + file_name[i] + ".png\'; plot \'" + export[i] + "\' using 1:2") for i in range(len(file_name))]

#selected = file_name #["serial_stream"] #file_name #["parallel_db_edt_breakdown_stream"]
#[subprocess.call(j.split("  ")) for i, j in plot_cmd if i in selected]
#print [j.split("  ") for i, j in plot_cmd if i in selected]

#color0 = "lt rgb \"red\""
#color1 = "lt rgb \"blue\""
#color2 = "lt rgb \"green\""
#color = [color0, color1, color2]
#prefix = plot_prefix + plot_dest + "all.png\';" + title + x_label + y_label 
#for j in [" plot \'" + export[i] + "\' " + color[i] + " using 1:2, " for i in range(len(file_name))]:
	#prefix += color[i]
#	prefix += j
#prefix[-1] = ';'
#all = plot_prefix + plot_dest + "all.png\';" + title + x_label + y_label + " plot \'" + "./results/ss.dat \' lt rgb \"blue\" using 1:2,  plot \'./results/ps.dat \' lt rgb \"green\" using 1:2;"
#print prefix
print all
#subprocess.call(["gnuplot", "-e", all])
# Need to take care of split_size

#Need to plot Native vs serial vs parallel vs parallel_db_edt_breakdown
