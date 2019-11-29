# RMRF Code

This program implements the RMRF algorithm and all results from the "Repeated Measures Random Forest (RMRF): Identifying Factors Associated with Nocturnal Hypoglycemia" manuscript from P. Calhoun, R.A. Levine, and J. Fan in Biometrics.

The RMRF algorithm and other various random forest algorithms are located in the "Functions" folder.

The results in Section 3 are located in the "Simulations" folder; the results and dataset in Section 4 are located in "Diabetes Study" folder.  The code does not need to be modified, although the code relies on the folder structure.  The working directory must be the path of the file opened.  The RMRF functions depend on several various R packages and loaded in the beginning of the file; tabulating results often imports other R packages too.  The code to generate some of the simulation results and construct the RMRF for the diabetes study takes a long time to run.  To reduce the computation time, the number of trees is reduced as described in the code comments.  The results in the manuscript are also saved and reproducible in the code.
