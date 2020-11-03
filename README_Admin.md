# bioinfo-bootcamp-fall-2020
Repo for bioinformatics bootcamp fall 2020

## Setting up the lab machines

Here's the protocol for setting up the lab machines:

1. IMS sets up an Azure lab template
2. Change the password to something easier to remember...
3. Start the template and connect via ssh
4. Install conda (saying 'yes' when prompted):

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

5. Use conda to install `salmon`, `refgenie`, and `sra-tools`:

```
conda install -c bioconda -c conda-forge salmon sra-tools refgenie
```

6. Pull the hg38 Salmon index

```
export REFGENIE='genome_config.yaml'
refgenie init -c $REFGENIE
refgenie pull hg38/salmon_partial_sa_index
```


7. Download SRA file and convert to fastq
```
prefetch SRR1039523
fasterq-dump -e 8 -p -O SRR1039523/ "/home/labadmin/SRR1039523/SRR1039523.sra"
```

8. Align to transcriptome using `salmon`
```
salmon quant -l A -1 "/home/labadmin/SRR1039523/SRR1039523.sra_1.fastq" -2 "/home/labadmin/SRR1039523/SRR1039523.sra_2.fastq" -i ~/hg38/salmon_partial_sa_index/default/ -p 8 -o Salmon.out/
```

9. To turn on all the machines before lecture (in powershell):

```
az account set --subscription research
az vm start --ids $(az vm list --resource-group m2600-workshop-basic-rg --query "[].id" -o tsv)
```

10. To export the mahine IP addresses to a csv:

```
$results = (az network public-ip list --resource-group m2600-workshop-basic-rg --query "[].{name: name, address: ipAddress}" -o json) |ConvertFrom-Json
$results|Export-Csv -Path vmlist.csv -NoTypeInformation
```

11. Sign in:

Username: Labadmin

Password: !L1nuxL@bT3$t!

12. To power off the machines

```
az vm deallocate --ids $(az vm list --resource-group m2600-workshop-basic-rg --query "[].id" -o tsv)
```
