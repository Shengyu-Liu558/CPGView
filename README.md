#                                         CPGView (Website: http://www.1kmpg.cn/cpgview/)

The goal of the entire project is to achieve genetic visualization, The entire program needs to be executed in the Linux environment. Below are the instructions for use of  CPGView:

### 1. Environment setting: 

​    1.1  The User need to download and install Perl under Linux

​			The following methods are for reference only.

```shell
wget http://www.cpan.org/src/5.0/perl-5.26.1.tar.gz
tar zxvf perl-5.26.1.tar.gz
cd perl-5.26.1
./Configure -de
make
make test
make installshel
```

​    1.2  Python 3 environment (Python 3.6 and above)

​			The user can click the link (https://www.python.org/) to select the appropriate python version and download.

​	1.3  Some related package in python

​			The user can execute the following command in python environment.

```shell
pip install Bio
pip install PyPDF2
pip install fitz PyMuPDF
```

​	1.4  R language environment

​			The user can click the link (https://www.r-project.org/) to select the appropriate R version and download.

​	1.5  Some related packages for visualization in R

​			1) ggplot2, gggenes, ape, dplyr, circlize, magrittr packages

```R
install.packages("ggplot2")
install.packages("gggenes")
install.packages("ape")
install.packages("dplyr")
install.packages("circlize")
install.packages("magrittr")
```

​			2) coRdon and genbankr packages

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("coRdon")
```

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("genbankr")
```



### 2. The process of implementation 

​	2.1  Upload the GenBank file

​	The user can upload the GenBank file to the Linux system, The suffix of the GenBank file are ".gb" or ".gbf". 

​    2.2  Change the location

​    The user needs to change the path to the location of R and Python in the user's local computer in "plasdrawmap.pl" file and "cpgview_wrapper.pl" file. The specific replacement location has been noted in the document.

​	2.3  Execute program

   "plasdrawmap.pl" is a perl implementation file(in the  "Linux"  folder), containing the phase of three image generated.  "plasdrawmap.pl" is a simple script provided using streaming API. It takes four arguments: the plasdrawmap.pl file, the input GeneBank file(Relative or absolute path), file name/ID and the name of output folder(Relative or absolute path). Since three jobs need to be running, this process maybe take 20 -40 seconds. Below is an example:

```shell
format: perl plasdrawmap.pl example_name result_folder example.gb
example: perl plasdrawmap.pl Arabidopsis Arabidopsis_folder sequence.gb
```

2.4  Prompt information

​			During the execution of "plasdrawmap.pl", some information will appear for users' reference, such as:

```
---  plasdrawmap completed! ---
```



### 3. The output folder

​			The output folder contains some images and some data folder. Users only need to pay attention to the files ending with "cc02.pdf", "cis.pdf" and "trans.pdf". Other files are generated in the middle for reference only.

### 4. Contact

​			Shengyu Liu - shengyuliu558@163.com 
​			Chang Liu - cliu6688@yahoo.com


