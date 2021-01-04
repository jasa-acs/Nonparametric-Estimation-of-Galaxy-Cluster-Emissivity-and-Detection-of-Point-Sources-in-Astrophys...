%First follow instructions in 'readme.txt' file

WavePath

%SIMULATION FOR TABLE 1 AND FIGURE 2
%Coment this loop to perform just the analysis%
for simoptions=1:24
	runSimTable1parallel
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analyseSimTable1
plotFigure2

%SIMULATION FOR TABLE 2 AND FIGURE 3
%Coment this loop to perform just the analysis%
for simoptions=1:16
	runSimFigure3parallel
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analyseSimFigure3
figure
plotFigure3 

%To produce latex version of tables 1 and 2, 
%use R software, and createTable1_2.R


%CODE FOR FIGURE 4 AND FIGURE 5 (REAL DATA)
%Coment this loop to perform just the analysis%
for simoptions=1:90
	runRealDataParallel
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analyseRealData
figure
plotFigure4
figure
plotFigure5


