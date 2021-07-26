## Table of contents  
1. [About this Repository](#About-this-Repository)
2. [Instructions to use the Jupyterhub on the VSC](#Instructions-to-use-the-Jupyterhub-on-the-VSC)
3. [scRNAseq start notebooks](#scRNAseq-start-notebooks)

## About this Repository  
This repository gives detailed instructions on how to use Jupyterhub on the [Vienna scientific cluster](https://vsc.ac.at/home/) with main focus on R applications. There are also two jupyternotebooks to get you started with your scRNAseq analysis.  
<br/><br/>
## Instructions to use the Jupyterhub on the VSC  
__step 1:__  if you affliated with the [Department of Dermatology](https://www.meduniwien.ac.at/web/en/about-us/organisation/university-departments/department-of-dermatology/?L=3) of the Medical University of Vienna and you want to use Jupyterhub via the VSC please do contact me: matthias.wielscher@meduniwien.ac.at. I will need your username and a cellphone number preferentially registered at an Austrian moblie network provider.  
__step 2:__ I will create a user on the VSC for you and we will have a brief meeting/call with a technical introduction and I will give an overview of code of conduct for using the VSC resources. Additionally, I strongly recomend attending [introductory classes](https://vsc.ac.at/research/vsc-research-center/vsc-school-seminar/) provided by the VSC team.   
__step 3:__ To start a Jupyternotebook go to here [https://vsc.ac.at/jupyterhub/hub/spawn](https://vsc.ac.at/jupyterhub/hub/spawn). You will be asked for your username a password, which will be provided by me. Then you will be asked to perform a 2-factor authentification with the cell phone number you provided in step1.  
__step 4:__ select VSC singularity image, then pick R notebook incase you want to use R. Basic/Data scienece/etc. notebooks are for python programming. You should look at screen that looks similar to that:  

<p align="center">
  <img src="pics/Screenshot_start_page.png" width="650" alt="accessibility text">
</p>
__step 4:__ Select your server specifications: run time, CPU, RAM memory


## scRNAseq start notebooks  
These are Jupyternotebooks for scRNAseq analysis. Comments and more detailed explainations are in the notebooks.  
>* this [notebook](singleSAMPLE_scRNAseq.ipynb) takes output from [10xcell ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview) and performs an initial exploratory analysis of the reads incl. doublet score calculation. Then it runs the Seurat SCT pipeline. The notebook als contains scirpt cells for differential gene expresssion analysis. It compares each cell population (cluster) against each remaining cluster (cell population) in the analysis. Resulting lists can be exported or can be further analysed with enrichR to gather information about possible cell type identities. Finally, cluster specfic markers are TIERed so you can see the top markers associtated to each cluster.
>* this [notebook](integrate_multiple10x_samples.ipynb) merges and integrates several 10x scRNAseq datasests together and then performs the same differential gene expression analysis as outlined above.  
