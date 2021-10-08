# Quick start

This is a short guide on how to set up your mavetools development environment. For this you will want your mavedb
development environment set up as well. Look [here](https://github.com/VariantEffect/mavedb/blob/main/DEVELOPERS.md) for 
instructions on how to do that. This guide assumes you will be using PyCharm. To start off, clone the repository and 
open it in PyCharm. This example will assume you name the cloned repository mavetools. PyCharm can also clone 
repositories for you if you supply your login credentials/access token. PLEASE make sure you switch to your feature 
branch.

Now, create a virtual environment in the root directory of the project:

```shell
python3 -m venv ./venv
```

Activate the virtual environment:

```shell
source ./venv/bin/activate
```

# Requirements

The following software is required:

- Python
- Git

Please install the following packages within your virtual environment:

- jupyter
- cfgv
- fqfa
- mavehgvs
- nbshpinx
- pandas
- renku-sphinx-theme
- Sphinx
- attrs

# Using Jupyter notebook

You can make requests via the API, Jupyter notebook examples can be found at 
`~/PycharmProjects/mavetools/docs/source` or on
[github](https://github.com/VariantEffect/mavetools/tree/two_tools/docs/source).

To run your Jupyter notebook from Pycharm, navigate to one of the .ipynb files in the directory described above and click 
the `Run All` (overlapped double green arrow) icon in the upper left of the editing window. You will notice  
a new tab called `Jupyter` available near your `Terminal` tab in Pycharm. Open this tab and select 
the `Server: mavetools` window. As described in the red text, copy and paste the `localhost` link into your browser. 
From here, you can follow the instructions within the notebook.

Note that, before running these scripts, you must have a local instance of MaveDB running.  