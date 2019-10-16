# 3.    Overview of Input

Walkabout requires information about the computational mesh, the flow solution, particle source locations, dispersion tensors, and numerical control parameters. Similar to the FEHM code, input is taken from a set of input files, the names of which are provided in the name file.

Conventions used here are as follows. Input parameter names, keyword types, and file types are shown in **bold font**. Literal character strings representing allowed values for keywords or filenames are shown in *italics*. The `monospaced font` is used when specifying input file formats. A line consisting of a single colon in the input block means that lines are skipped. An ellipsis (...) indicates that the item is to be repeated. Any text following an exclamation point in an input block is to be regarded as a comment or explanation.

## 3.1  Running Walkabout

Walkabout is started by typing the name of the executable file from a command line interface – no graphical user interface is provided. There is one optional command line argument giving the relative path to the **name** file.

## 3.2  The Name File

Names for the various input and output files are contained in the **name** file. The **name** file itself may be specified on the command line when invoking Walkabout. Otherwise it defaults to *walkabout.files* in the current working directory.

The format for the **name** file is as follows

**`filetype:filename`**<br/>
`: ! repeat as many times as needed`

The colon separating the **filetype** keyword and the **filename** is required. The **filetype** keywords may come in any order. A line that does not start with a valid **filetype** is ignored, thus allowing the user to add as many comments as desired to the file. Filenames may be given as full or relative filepaths as defined by the operating system. Valid **filetypes** are as follows.

*control* – specifies the control file, which contains control parameters and the dispersion tensor information. Optional. Defaults to *control.dat*

*fehmn* – specifies the FEHM geometry file. Required.

*stor* – specifies the FEHM stor file. Required.

*ealist* – specifies the element adjacency list. Required.

*fin* – specifies the FEHM fin file, which contains liquid fluxes. Required.

*avs* – specifies an FEHM AVS file, which must contain liquid saturation, porosity and liquid density. Required.

*cbound* – specifies a file containing a list of nodes that lie on boundaries that are closed to transport file. Optional.

*trajout* – specifies the output filename for trajectory output. Optional. Defaults to *walkabout.out*.

*sptr2* – specifies the output filename for use as PLUMECALC input. Optional. Defaults to *walkabout.sptr2*.

Details for each file type are given in the appendix.