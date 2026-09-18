import os, sys

#description:
#python function writing string to file

readfile = sys.argv[1]
output_dir = sys.argv[2]
string = sys.argv[3]
initiliaze = sys.argv[4]

if not os.path.isdir(output_dir): os.mkdir(output_dir)

mode = "a+"
if initiliaze == "true": mode = "w+"
fileread = open(output_dir+"/"+readfile,mode)
fileread.write(string)

fileread.close()