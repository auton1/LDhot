LDhot
=====

Detect recombination hotspots using population genetic data.

#Installation

After downloading, switch to the download folder and type:
```
make
```

If you're lucky, this will compile without errors. However, note you may need a compiler that supports the C++11 standard. If you see lots of errors, you may need to upgrade your compiler.

On some systems, you can compile with multi-threading turned on, which results in a signficant reduction in runtime. To do this, type:
```
make MULTI=1
```


#Usage

Two programs are provided. They are called as follows.

```
./ldhot --seq <seq_file> --loc <loc_file> --lk <lk_file> --res <res_file> --nsim 1000 --out <output_prefix>
./ldhot_summary --res <res_file> --hot <hotspot_file> --out <output_prefix>
```

The seq\_file, loc\_file, lk\_file, and res\_file are all derived from LDhat (http://ldhat.sourceforge.net/). The --nsim parameter controls the number of simulations used within the method, with at least 1000 being recommended.

The ldhot program produces an output file of the form \<output\_prefix\>.hotspots.txt. This file contains the details of the windows tested for the presence of a hotspot. This file can be treated as the final product, or further summarized using the ldhot\_summary program. This program uses the output of ldhot in conjunction with the recombination rate estimates to find hotspot boundaries. The output of this program can be found in \<output\_prefix\>.hot\_summary.txt.
