# eFindSite
eFindSite is a tool that predicts binding pockets, residues and ligands from a given protein structure by threading methods. 

This README file is written by Snigdha Thumma.

If you find this tool useful, please cite these following papers:

1.	Brylinski M, Feinstein WP. (2013) eFindSite: Improved prediction of ligand binding sites in protein models using meta-threading, machine learning and auxiliary ligands. *J Comput Aided Mol Des* **27** (6): 551-67.

2.	Feinstein WP, Brylinski M. (2014) eFindSite: Enhanced fingerprint-based virtual screening against predicted ligand binding sites in protein models. *Mol Inf* **33** (2): 135-50.  

# Prerequisites:
- GCC compiler version 4.8.5+
- Perl version 5.16.3+


------

# Example:
1. First you have to download eFindSite from here `https://github.com/michal-brylinski/efindsite` and unzip: 
     - `[local]$ unzip efindsite-master.zip`
     - `[local]$ cd efindsite-master/`
     - `[efindsite-master]$ make`
     - `[efindsite-master]$ ls bin/`
     - `[bin]$ ./efindsite`


2. Then install the eFindSite template libraries from here `https://osf.io/mp343/`*( download the latest available library )* and unpack the tarball:
     - `[local]$ mkdir libraries`
     - `[local]$ tar xzf efindsite-lib-2018-04.tar.gz -C /usr/local/libraries/`
     - `[local]$ tar xzf efindsite-mod-2016-06.tar.gz -C /usr/local/libraries/`

3. Before running eFindSite, you need to set several environmental variables. 
     - eFindSite
       - `[home]$ export EF_LIB=/usr/local/libraries/efindsite-lib-2018-04`
       - `[home]$ export EF_MOD=/usr/local/libraries/efindsite-mod-2016-06`

4. After the environmental variables are set, it’s time to run the eFindSite and efindsite_screen
     - `[library]$ /usr/local/efindsite-master/bin/efindsite -s 13gsA.pdb -t 13gsA-efindsite.lst -e 13gsA.profile -o 13gsA-efindsite_test`
     - `[library]$ /usr/local/efindsite-master/bin/efindsite_screen -p 13gsA-efindsite.pockets.dat -s escreen-keggcomp-mar2012.gz -o 13gsA-escreen-keggcomptest`
     - `[library]$ /usr/local/efindsite-master/bin/efindsite_screen -p 13gsA-efindsite.pockets.dat -s escreen-zinc12_nr-mar2012.gz -o 13gsA-escreen-zinc12_nrtest`

------

Mandatory arguments for eFindSite

|Parameter  |  Optional   |  Name |  Description   |
|:---:|:---:|:---:|:---:|
|    -s            |         N          |   input_file           | Is the target protein structure in PDB format; you can use either experimental structure or protein model |
|     -t            |          N         |     template-fun  | Is a text file that contains information on protein templated identified for your target e.g. by eThread |
|     -i            |          N         |     sec_struct.ss  | Is a secondary structure profile constructed for the target by PSIPRED |
|     -e           |          N         |     seq.prf            | Is a sequence profile of the target; it can be generated by PROFILPRO |
|     -o           |          N         |     output_name  | Is used to save five output files with different extensions that store dinging site prediction results |

Mandatory arguments for efindsite_scr

|Parameter  |  Optional   |  Name |  Description   |
|:---:|:---:|:---:|:---:|
|     -p          |      N   |    pocket.dat    |     pocket.dat is a file containing the predicted binding pockets by eFindSite  |
|     -s          |      N   |   cmp_lib          |      cmp_lib is a compound library used for virtual screening           |
|     -o          |    N     |    output_file    |   output_file will contain ranked library compounds assigned with a Tanimoto score and Z-score  |

------

- If eFindSite runs successfully, it should yield 5 output files with user defined input name.
** For example: the user defined input file, 13gsA-efindsite_test, has 5 outputs as seen below*
     - 13gsA-efindsite_test.pockets.dat (detailed info on predicted pockets)
     - 13gsA-efindsite_test.pockets.pdb (predicted pockets in PDB format)
     - 13gsA-efindsite_test.alignments.dat (structure alignments of templates into the target in FASTA format)
     - 13gsA-efindsite_test.templates.pdb (template structures aligned onto the target in PDB format)
     - 13gsA-efindsite_test.ligands.sdf (extracted binding ligands in SDF format)

- Open the output files and check to make sure the files contain data.
- If the output files are empty, go back and check the script.
