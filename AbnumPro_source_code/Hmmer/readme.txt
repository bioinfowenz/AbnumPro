
			      HMMER 3.0 for Windows

This suite of HMMER 3.0 modules has been compiled under CYGWIN environment 
and tested to run properly under Windows XP and Windows 7 operating systems.

In order to run HMMER programs, one need to place two dll files, 
cygwin1.dll and cyggcc_s-1.dll into directory, where those files can be seen by Windows, i.e.
place it in windows/system32 folder, or append the path to those files to PATH 
environment variable.

Note, that when you run HMMER modules and indicate file path as a parameter, you should use path either 
relative to the place, where HMMER modules located or use notation /cygdrive/c/mypath/myfile instead of
c:\mypath\myfile. Otherwise, you will every time get an annoing warning like this below:
 
cygwin warning:
  MS-DOS style path detected: c:\temp\file
  Preferred POSIX equivalent is: /cygdrive/c/temp/file
  CYGWIN environment variable option "nodosfilewarning" turns off this warning.



Sergey Smirnov
NCBI/NLM/NIH

email: sergey@ssmirnov.com; smirnov@ncbi.nlm.nih.gov  