LDhot
=====

A program to detect recombination hotspots using population genetic data.

##Installation

After downloading, switch to the download folder and type:
```
make
```

If you're lucky, this will compile without errors. However, note you may need a compiler that supports the C++11 standard. If you see lots of errors, you may need to upgrade your compiler.

On some systems, you can compile with multi-threading turned on, which can result in a signficant reduction in runtime. To do this, type:
```
make MULTI=1
```


##Basic Usage

Two programs are provided. The main program is called as follows.

```
./ldhot --seq <seq_file> --loc <loc_file> --lk <lk_file> --res <res_file> --nsim 1000 --out <out_prefix>
```

The **seq\_file**, **loc\_file**, **lk\_file**, and **res\_file** are all derived from [**LDhat**](http://ldhat.sourceforge.net/). The --nsim parameter controls the number of simulations used within the method, with at least 1000 simulations being recommended. A complete option list is given below. 

The **ldhot** program produces an output file of the form **\<output\_prefix\>.hotspots.txt**, which contains the details of the windows tested for the presence of a hotspot. 
This file can be treated as the final output, or further summarized using the **ldhot\_summary** program. This is a simple program which combines windows called as significant by the main **ldhot** program. It is called as follows.

```
./ldhot_summary --res <res_file> --hot <hotspot_file> --out <out_prefix>
```

The output of this program can be found in **\<output\_prefix\>.hot\_summary.txt**.

A more complete example of the usage of LDhat and LDhot, with both input and output files, can be found in the **example** folder.

##Option List

###ldhot

The **ldhot** program takes the following parameters.

####Required Parameters:
* --seq <filename> : Input LDhat-format sequence file.
* --loc <filename> : Input LDhat-format positions file.
* --lk <filename>  : Input LDhat-format likelihood lookup file.
* --res <filename> : Input recombination rate estimates in same format as LDhat 'stat' output.

####Important Parameters:
* --out <prefix>   : Prefix for output files (default: out).
* --nsim <int>     : Maximum number of simulations to use (default: 100 but at least 1000 recommended).

####Other Parameters:
* --startpos <double>   : Start position in kb.
* --endpos <double>     : End position in kb.
* --step <double>       : Step size (in kb) between tested windows (default: 1).
* --windist <double>    : Define background window as +/- *windist* kb of hotspot center (default: 50).
* --hotdist <double>    : Define hotspot window as +/- *hotdist* kb hotspot center (default: 1.5).
* --seed <int>          : Random seed.
* --nofreqcond          : Turn off frequency conditioning.
* --lk-SNP-window <int> : Number of SNPs over which to calculate the composite likelihood (default: 50).

###ldhot_summary

The **ldhot_summary** program takes the following parameters.

####Required Parameters:
* --res <filename> : Input recombination rate estimates in same format as LDhat 'stat' output.
* --hot <filename> : Input hotspot file from LDhot.

####Other Parameters:
* --out <prefix>   : Prefix for output files (default: out).
* --sig <double>      : Significance cutoff for calling a hotspot (default: 0.001).
* --sigjoin <double>  : Significance cutoff for merging hotspot windows (default: 0.01).