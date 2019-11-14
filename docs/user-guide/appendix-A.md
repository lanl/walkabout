# Appendix A

## NAME File

Appendix A describes the Walkabout input and output files. The names for these files are contained in the Walkabout name file. 


The **name** file may be specified on the command line when invoking Walkabout. Otherwise it defaults to **walkabout.files** in the current working directory.


The format for the **name** file is as follows:
```
filetype:filename
: ! repeat as many times as needed
```
A line that does not start with a valid filetype is ignored, thus allowing the user to add as many comments as desired to the file. Filenames may be given as full or relative filepaths as defined by the operating system. Each file is described in the section indicated.

*Note, some filetypes have default names as indicated below, and therefore are optional for the name file.*


### Required filetypes for NAME file:

* **`control`** specifies the control file, which contains control parameters and the dispersion tensor information.  Default name is control.dat [`Appendix A.1 Control`](appendix-A1.md)

* **`fehmn`**  specifies the geometry in FEHM mesh format. [`Appendix A.2 Geometry`](appendix-A2.md) 

* **`stor`** specifies the geometric coefficients in stor format. [`Appendix A.2 Geometry`](appendix-A2.md)

* **`ealist`** specifies the element adjacency list. [`Appendix A.2 Geometry`](appendix-A2.md)

* **`fin`** or **ama** specifies the FEHM fin file or AMANZI ama file which contains liquid fluxes. [`Appendix A.3 Model`](appendix-A3.md)

* **`avs`** specifies an AVS format file with node properties liquid saturation, porosity and liquid density. [`Appendix A.3 Model`](appendix-A3.md)


### Optional filetypes for NAME file:

* **`cbound`** this file may be used to list nodes that lie on outside boundaries that are closed to transport. [`Appendix A.2 Geometry`](appendix-A2.md)

* **`trajout`** specifies the output filename for trajectory output.  Default name will be walkabout.out [`Appendix A.4 Output`](appendix-A4.md)

* **`sptr2`** specifies the output filename to use for the PLUMECALC file.  Default name will be walkabout.sptr2 [`Appendix A.4 Output`](appendix-A4.md)


The following sections Appendix A.1 to A.4 describe these Walkabout filetypes and their formats.


