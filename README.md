# csv_notebook
Python notebooks to read RD meas from yield_ratio*.json file and calculate CSV from the json file
- Yield_ratio_check.ipynb is the one compare RDmeas as a function of z for one Q2 and xbj. It calculated CSV from RDmeas(data) and some existing model, cteq6l1 for pdf, fdss for ff

- Yield_ratio_zx_check.ipynb only compare RDmeas as a function of x for one Q2 and z

- Yield_ratio.ipynb tried to fit all the rdmeas without ff model to get our own CSV and FF

- CSV_calculation_Q2center.ipynb: calculate csv from existing FF for Three Q2 center, including FFs from my fitting from yield_ratio_morepoints.ipynb, FFs params are saved in CSVs_2dfit.json with keyword "7ass". calculated CSV are also saved in CSVs_2dfit.json, which is used to plot along the fitting from yield_ratio_morepoints.ipynb
- CSV_calculation_Q2center_rho.ipynb: calculate csv from existing FF for three Q2 center, including FFs from mt fitting from yield_ratio_morepoints_rho.ipynb, FFs params are saved in CSVs_2dfit.json with key "7ass_rho". calculated CSv are also save in CSvs_2dfit.json, which is used to plot along the fitting from yield_ratio_morepoints_rho.ipynb

-CSV_calculation_Q2cut.ipynb: calculate csv from existing FF with Three hard Q2 cut. Since we are using Q2 corrected, this one is not necessary, abandon 

-Yield_ratio_morepoints.ipynb: tried to fit all the rdmeas without ff model to get my own CSV and FF, the difference with yield_ratio.ipynb is that, this one used the json input with more points, same kinematic point with different rungroup not combined
-Yield_ratio_morepoints_pandas.ipynb: same as previous one, but using csv.csv as input instead of json file. So it's easier to save residual information. The fitting results is slightly different, since the json input has more precision than csv.csv input
-Yield_ratio_morepoints_rho.ipynb: tried to fit all the rdmeas after rho subtraction without ff model to get my own CSV and FF
-Yield_ratio_morepoints_Arho.ipynb: tried to fit all the rdmeas without ff model to get my own CSV and FF, here the rho subtraction is down with one param a_rho
-Yield_ratio_morepoints_Arho_pandas.ipynb: tried to fit all the rdmeas without ff model to get my own CSV and FF, here the rho subtraction is another input param in least square. The input data is from csv.csv

-Yield_ratio_morepoints_scipy.ipynb, Yield_ratio_scipy.ipynb: instead of using iminuit to minimize, I used scipy. Not much difference, not updating for a while

#Feb11,2022
-H2runs.ipynb, read H2 and D2 data from csv_H2.csv, plot diff ratio and sum ratio

#I uploaded several different results_W2Wp2_#_# saves all json and csv file with different W2 and Wp2 cut , by soft link 
ln -s results_W2Wp2_#_# results
the above ipynb, yield_ratio_morepoints_pandas.ipynb, CSV_calculation_Q2center.ipynb and H2_runs.ipynb will read the json and csv file from results, and print the result to the results
