# Getting started

## Quick start

There are several ways you can work with SWAP. You can use the traditional way, where you download the model exacutable from the distributing website, or you can choose to use one of the interfaces (currently available in R and Python).

Follow one of the following ways to get started with SWAP :rocket:

### :material-download: _Manual download_

Go to the [Alterra website](https://www.swap.alterra.nl/). Under Sownloads tab, you will find a .zip file containing the model executable and some additional resources.

Once you download the file, unzip it at a desired location. To start with, you can navigate to the directory titled 'cases'. Inside, you will find several examples of SWAP usage. Navigate to one that seems the most interesting to you. You will find a `run_swap.cmd` file there. Double click on it and wait for the programme completion.

You will notice that a bunch of new files appeared. These are the model outputs. They are simple ascii (or CSV) files which you can subsequently load into a data analysis software of your choice.

### :simple-github: _GitHub_

!!! note

    This is an option that we should perhaps include. It would involve cloning the code and compiling it on users' machines. This would also support forking and community driven feature development.

### :simple-r: _packages_

You can consult the [documentation of rSWAP](https://moritzshore.github.io/rswap/) to see how you can install and run SWAP models from RStudio.

If you follow the steps of The Traditional Way, you will also find another R package called SWAPTools, distributed along with the SWAP model executable in the `rsoftware/` directory. In the `docs/` directory, you will find instructions on how to install it in RStudio. You can decide which package fits better your needs.

### :material-language-python: _Python package_

If you are a Python programmer, you can use pySWAP: Python wrapper for SWAP model. You can follow [the documentation](https://zawadzkim.github.io/pySWAP/) to install the package and run the first model.
