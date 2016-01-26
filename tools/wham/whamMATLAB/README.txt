These scripts will post process the data from Aaron Wenger's RE program.

First:
Copy the entire wham folder into the directory where all of the RE output files are.

Second:
Run the makeDihedralEnergyFiles.sh script to extract the necessary datafiles.

Third:
IN MATLAB run the plot_re.m script which will call the remaining MATLAB scripts.
(WHAM.m will be run and the plot generated)
	example: matlab -nodisplay < plot_re.m > matlabscreen.out &
