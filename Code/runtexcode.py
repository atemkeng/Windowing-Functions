import numpy as np
import pickle
def run(filename):
	f = open(filename, "r")
	Hires = []; Lores = []; Noise = []; Window = [];
	try:
		while True:
			list =  pickle.load(f);
			Hires.append(list.get("HiresMS"));
			Lores.append(list.get("LowresMS"));
			Noise.append(list.get("Noise"));
			Window.append(list.get("window"))
	except:
		print "\n"
	n = len(Window)
	strgg = ""
	for i in range(n): 
		aa = r"%s&Nx(%ix%i)&%.3f&%.3f&-&-&-\\"\
		%(Window[i],Lores[i][0],Lores[i][1],Noise[i]*10**3,Noise[i]/Noise[-1])
		strgg = strgg+""+aa
	cccc =r"""
	\hspace{-0.0cm}\begin{tabular}{|l||l|l|l||l|l|l|}
	\hline
	 \textbf{Method} &\multicolumn{3}{l|}{\textbf{Standart Technique}}&\multicolumn{3}{l|}{\textbf{Baseline dependent Technique}}\\
	\cline{2-7}
	 &comp. rate&rms(mJy)&rms ration&comp. rate&rms(mJy)&rms ration\\
	\hline\hline
	%s
	\hline
	Hires dataset  &Nx(%sx%s)\\
	\cline{1-2}
	\end{tabular}
	"""%(strgg,str(Hires[0][0]),str(Hires[0][1]))
	return cccc

if __name__ == "__main__":
	import os;
	import sys;
	arg = sys.argv[1]
	ss = run(arg); 
	latex = r"""
\documentclass[a4paper]{book}
\usepackage{dcolumn}
\author{M. Atemkeng}
\begin{document}
\section*{Compression rate and Noise}
%s
\end{document}"""%ss
	#np.save("test1.tex",latex)
	f = open('%s.tex'%(arg[0:-5]), 'w')
	f.write(latex)
	f.close()
	command = "pdflatex -shell-escape %s.tex"%(arg[0:-5])
	os.system(command)




