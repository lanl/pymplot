
# `pymplot`: A lightweight matplotlib-based application package for plotting 1D, 2D and 3D scalar datasets

## Description

**pymplot** is a lightweight Python application package for efficiently generating high-quality plots for multi-dimensional scalar data. See the paper draft in doc for more details. 

LANL C number: C20105

Author: Kai Gao, <kaigao@lanl.gov>

## Installation

	cd src
	ruby install.rb
	
The installed links to executables are installed at $HOME/bin, with names `x_showmatrix`, `x_showcontour`, `x_showwiggle`, etc. Please make sure that $HOME/bin is in your system path. 

You also need to have `numpy`, `scipy`, `matplotlib`, and `pyvista` in your Python distribution. 

Currently, Windows platforms are not supported due to technical reasons. 

## Documentation
You can type `x_show...` in the commond line to check the meaning of arguments. 

## Examples

Please refer to the ruby script (`example/test.rb`) to reproduce some of the examples in the paper draft.

## License

**pymplot** is distributed under the `BSD license`. See details in LICENSE. 


## Reporting Bugs or Suggest Improvements

Report bugs [here](https://github.com/lanl/pymplot/issues). If you are reporting a bug, please include:

* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

The author welcomes functionality improvement/extension suggestions. If you are suggesting an improvement/extension, please include:

* A detailed description of the suggested functionality(ies).

