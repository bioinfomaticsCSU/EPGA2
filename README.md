License
=========

Copyright (C) 2014 Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


De Novo Assembler
=================

EPGA2 updates some modules in EPGA which can improve memory efficiency in genome asssembly.
The read library for EPGA2 should be paired-end reads. Read length shorter than 50bp and coverage larger than 100.

1)Installing.

EPGA is written C++ and therefore will require a machine with GNU C++ pre-installed.

Create a main directory (eg:EPGA2). Copy all source code to this directory.

2)Running.

Run command line: 

ulimit -n 1100 //this command is used for BCALM 

g++ main.cpp -o EPGA -lpthread 

perl EPGA.pl library.txt kmerLength threadNumber

Information about read libraries are stored in the file library.txt.
Each line represents one read library.
The first column is the first mate read file (*.fastq), the sencond column is the second mate read file (*.fastq), the third column is insert size of read library, the fourth column is standard deviation of insert size, the fifth column represents whether the read library is mate-paired (0 denotes paired-end reads, 1 denotes mate-paired reads).

kmerLength is one integer (<32) shorter than read length which is used for building De Bruijn graph.

threadNumber is thread number of program.

3)Output.

There are two files "contigSetLong.fa" and "scaffoldLong.fa" corresponding to final contigs and scaffolds.

