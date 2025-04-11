# Data

## Small molecule-miRNA pairs

There are 2 small molecule-miRNA (SM-miR) pairs as follows.

| No. | small molecule |    mirna     |
|:---:|:--------------:|:------------:|
|  1  |      5743      | MIMAT0000539 |
|  2  |   84093        | MIMAT0009203 |
|  2  |     148124     | MIMAT0000098 |
|  2  |      5757      | MIMAT0000066 |

## miRNAs

:::{note} Fasta
:class: dropdown
Protein sequences in [the Fasta format](https://en.wikipedia.org/wiki/FASTA_format) are required. The file extension must be `.fasta` for recognition of the software.
:::

## Small molecules

Drutai reads the smile strings of small molecules as follows. This needs to be separated by `tab`.

| No. | small molecule |            smile             |
|:---:|:--------------:|:----------------------------:|
|  1  |    5743        |  C[C@@H]1...=O)CO)O)C)O)F)C  |
|  2  |     84093      | CC1=CN(C(=...].[NH2-].[Pt+2] |
|  2  |     148124     | CC1=C2[C@H...OC(=O)C)O)C)O   |
|  2  |      5757      | C[C@]12CC[C@...=C3C=CC(=C4)O |


## Example data

Users can download our example data this way.

::::{tab-set}
:::{tab-item} Code
:sync: tab1
```{code} python
import deepsmirud

deepsmirud.predict.download_data(
    url='https://github.com/2003100127/deepsmirud/releases/download/example-data/example_data.zip',
    sv_fpn='../../data/deepsmirud/example_data.zip',
)
```
:::
:::{tab-item} Output
:sync: tab2
```{code} shell
 ____                                _      _   _ ____  
|  _ \  ___  ___ _ __  ___ _ __ ___ (_)_ __| | | |  _ \ 
| | | |/ _ \/ _ \ '_ \/ __| '_ ` _ \| | '__| | | | | | |
| |_| |  __/  __/ |_) \__ \ | | | | | | |  | |_| | |_| |
|____/ \___|\___| .__/|___/_| |_| |_|_|_|   \___/|____/ 
                |_|                                     

05/04/2025 12:00:17 logger: =>Downloading starts...
05/04/2025 12:00:18 logger: =>downloaded.
```
:::
::::

