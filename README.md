R scripts supporting the manuscript "The Effect of Screening for Publication Bias on the Outcomes of Meta-Analyses" by Haben Michael (hmichael@math.umass.edu).

``sim.R`` implements the simulation described in the simulation section of the manuscript. 

``process.R`` generates the tables in the manuscript from the output of ``sim.R``.

``utils.R`` contains helper functions called by ``sim.R``.

``submit.sh`` is a shell script submitting ``sim.R`` to a job scheduler (slurm).

``main.R`` is an ```R``` script to generate the tables using a single thread.


To replicate the simulation in Section 5 of the manuscript, first run
the script ```submit.sh``` from a shell with ```Rscript``` in the
path.  In ```submit.sh```, setting ```theta=0``` gives Table 1a and
```theta=0.2``` gives Ta- ble 1b. After execution completes the
working directory will contain ```R``` save files containing the results of
the simulation. Running ```process.R``` from an ```R``` session in the same
directory will then output the tables as ```Latex``` files.
