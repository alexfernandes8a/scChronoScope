# scChronoScope

An interactive Dash application that offers comprehensive visualization of gene expression dynamics in single-cell RNA sequencing data, capturing changes over time and across diverse cell types.

## Running

### Set Up the Appropriate Python Environment

- **Tested with Python version 3.13.2**

```shell
python3 -m venv ENV
source ENV/bin/activate
pip install -r requirements.txt
```

### Run the App
Make sure you have all the files in the correct folders and that you are in the correct environment. Then run:

```shell
python scChronoScope.py
```
Now, simply copy/paste/open the link in your preferred browser.



#### In case you need to set a specific port and binding host, run (optional):
```shell
python scChronoScope.py --host 127.0.0.1 --port 8050
```

That's it!




## OPTIONAL | Jupyter Lab Notebook Version:
If you prefer to run the application using Jupyter Lab, follow these steps in the same ENV:

1. Start Jupyter Lab:
```shell
jupyter lab
```

The Jupyter Lab interface should open in your browser. Navigate to and open the 'scChronoScope.ipynb' notebook.

Enjoy exploring the interactive visualizations!
