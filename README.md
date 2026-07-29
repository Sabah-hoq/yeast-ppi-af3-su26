# yeast
---
### Background 
Analyzing data from [GitHub repo](https://github.com/jurgjn/pooled-ppi-yeast). 

#### __Docker Data__
To pull data, run this line of code (approx. 30 GB). More information about the Docker Image can be found [here](https://hub.docker.com/r/jurgjn/pooled-ppi-yeast/tags). This has all the data needed to perform analysis 

```
docker run -p 8501:8501 jurgjn/pooled-ppi-yeast:v26.1
```
After running, the raw data of the proteins is most likely stored in `predictions-db`; the data I will extract/copy is `summary_confidences.parquet` and `summary_pairs.parquet` which are in `workspaces/data`.


_The lines below are for windows PowerShell_
```{shell}
> New-Item -Path "C:\yeast_data" -ItemType Directory -Force
> docker cp <container_name>:/workspace/data/summary_pairs.parquet C:\yeast_data\summary_pairs.parquet
> docker cp <container_name>:/workspace/data/summary_confidences.parquet C:\yeast_data\summary_confidences.parquet
> explorer C:\yeast_data
```
MacOS Code:
```{bash}
> mkdir -p ~/yeast_data
> docker cp <container_name>:/workspace/data/summary_pairs.parquet ~/yeast_data/summary_pairs.parquet
> docker cp <container_name>:/workspace/data/summary_confidences.parquet ~/yeast_data/summary_confidences.parquet
> open ~/yeast_data
```

For linux, to open the folder do:

```{bash}
xdg-open ~/yeast_data
```

* yeast_data: A folder I made locally; replace with the folder you want the data in
Basic analysis of data will be in: 

If running into issues where the protein interface won't load after running the Dockerfile, try:

```{shell}
docker run -d -p 8501:8501 --name <conatiner_name> jurgjn/pooled-ppi-yeast:v26.1
```
Comparing the STRING Scores with the AF3 ipTM scores is [here](https://github.com/Sabah-hoq/yeast/blob/main/notebooks/PR_ROC/graph.ipynb)

#### __Getting STRING Data__
As for right now, the `.gitignore` file doesn't track data files as they are large, but we want to make sure all of our data is stored in a data folder. 
To get STRING files run this in the scripts: [download STRING files](https://github.com/Sabah-hoq/yeast/blob/main/scripts/string__downloader.py)

The files used can found down below:
[Yeast](https://string-db.org/cgi/download?sessionId=boDI8ehQpWOh&species_text=Saccharomyces+cerevisiae)

The code from this file was copied from this [script](https://github.com/jurgjn/pooled-ppi/blob/main/src/pooled_ppi/string_db.py)

*Clicking on the links down below will automatically download the files*

* [protein alias](https://stringdb-downloads.org/download/protein.aliases.v12.0/4932.protein.aliases.v12.0.txt.gz)
* [protein information](https://stringdb-downloads.org/download/protein.info.v12.0/4932.protein.info.v12.0.txt.gz)
* [protein physical links detailed]([https://stringdb-downloads.org/download/protein.links.detailed.v12.0/4932.protein.links.detailed.v12.0.txt.gz](https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/4932.protein.physical.links.detailed.v12.0.txt.gz))
---
### Understanding ipTM & STRING scores

### Files 

**Cytoscape**

There is one file already in the `cytoscape` folder, `yeast_network.cys`. This file when downloaded, you can view the pathways of the yeast PPI from AF3. However, the datasets used have been trimmed/filtered, as the original unfiltered dataset had 7.6 million + edges, and Cytoscape does not have the memory capabilities for this. So, the script, [`cytoscape.py`](xdg-open ~/yeast_data) grabs a file from the data folder and and filters the ipTM and STRING scores where it is at a signigicant threshold. This is done for both the ipTM mean and best scores, and then visualzied using Cytoscape. 

