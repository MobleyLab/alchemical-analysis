An implementation of the recommended practices for analyzing alchemical free energy calculations, as described in Klimovich et al., JCAMD 29:397-411 (2015).

There are 11 subdirectories in `data/` each containing `ti00[23].out`.  The output files were created with the command in run.sh.

The argument that follows the `-a` flag should not necessarily be in all capitals; any combination of lower- and upper-case letters is OK.

`data` is the path to the directory with the data files.

The `-p` flag seeks for the prefix of the data file. If the data files are in multiple subdirectories,
the name of those subdirectories (in a form of the glob pattern) should preface the file prefix (like `ti*/ti` above).

The dataset contained in the `data/` directory was generated from the files obtained from a 11-window simulation of the electrostatic solution-phase methanol-to-methane transformation run by Hannes Loeffler at STFC, UK.
