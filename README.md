## Description
**pymplot** is a lightweight Python application package for efficiently generating high-quality plots for multi-dimensional scalar data.  

The work is under LANL open source approval reference C20105.

# Reference
Please refer to the [paper draft](doc/paper.pdf) for details. 

## Installation

	cd src
	ruby install.rb
	
The installed links to executables are installed at $HOME/bin, with names `x_showmatrix`, `x_showcontour`, `x_showwiggle`, etc. Please make sure that $HOME/bin is in your system path. 

You also need to have `numpy`, `scipy`, `matplotlib`, and `pyvista` in your Python distribution. 

## Documentation
You can type `x_show...` in the commond line to check the meaning of arguments. 

## Examples
Please refer to the [test.rb](example/test.rb) to reproduce some of the examples in the paper draft.

## License
`pymplot` is distributed under the `BSD license`. See details in LICENSE. 

# Author
Kai Gao, <kaigao@lanl.gov>

We welcome feedbacks, bug reports, and improvement ideas on `pymplot`. 

If you use this package in your research and find it useful, please cite it as

* Kai Gao, Lianjie Huang, 2021, Pymplot: An open-source, lightweight plotting package based on Python and matplotlib , url: [github.com/lanl/pymplot](github.com/lanl/pymplot)
