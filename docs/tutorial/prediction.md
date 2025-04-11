# Prediction

## Overview

You need to decompress the `example_data.zip` file in your preferred folder, e.g., `deepsmirud/`.

There are three required files.

* `br_sm_mirna.txt`
* `br_smile.txt`
* `path` to `*.fasta`

We display 4 pairs of small molecules and miRNAs to predict their possible regulation types.

```{code}
:linenos:
:emphasize-lines: 2,3,4,5
sm	mirna
5743	MIMAT0000539
84093	MIMAT0009203
148124	MIMAT0000098
5757	MIMAT0000066
```

:::{caution} Model declaration
* For TF `1.15x` models, including `alexnet`, `birnn`, `rnn`, and `seq2seq`, you need to declare them like `alexnet/alexnet`.
* For TF `2.14x` models, including `cnn`, `lstmcnn`, `dsconv`, `convmixer64`, `mobilenet`, `resnet_prea18`, `resnet_prea50`, and `seresnet`, you need to declare them like `lstmcnn`.
:::

## Python

We access **DeepsmirUD** by defining the following parameters. 

:::{seealso}
We have `method` defined as **AlexNet**, **BiRNN**, **RNN**, **Seq2Seq**, **CNN**, **ConvMixer64**, **DSConv**, **LSTMCNN**, **MobileNet**, **ResNet18**, **ResNet50**, or **SEResNet**. Please see also the DeepsmirUD paper (@Sun2023deepsmirud) for method details.
:::

```{code} python
params = {
    'br_fpn': '../../data/deepsmirud/example_data/br_sm_mirna.txt',
    'smile_fpn': '../../data/deepsmirud/example_data/br_smile.txt',
    'fasta_fp': '../../data/deepsmirud/example_data/',
    'method': 'LSTMCNN',
    'model_fp': '../../data/deepsmirud/model/lstmcnn',
    'sv_fpn': '../../data/deepsmirud/example_data/pred.deepsmirud',
}
```

::::{tab-set}
:::{tab-item} Command
:sync: tab1
```{code} python
import deepsmirud

deepsmirud.predict.sm_mir_regulation_type(
    br_fpn=params['br_fpn'],
    smile_fpn=params['smile_fpn'],
    fasta_fp=params['fasta_fp'],
    method=params['method'],
    model_fp=params['model_fp'],
    sv_fpn=params['sv_fpn'],
)
```
:::
:::{tab-item} Output log
:sync: tab2
```{code} shell
 ____                                _      _   _ ____  
|  _ \  ___  ___ _ __  ___ _ __ ___ (_)_ __| | | |  _ \ 
| | | |/ _ \/ _ \ '_ \/ __| '_ ` _ \| | '__| | | | | | |
| |_| |  __/  __/ |_) \__ \ | | | | | | |  | |_| | |_| |
|____/ \___|\___| .__/|___/_| |_| |_|_|_|   \___/|____/ 
                |_|                                     

05/04/2025 12:03:09 logger: =>Prediction starts...
05/04/2025 12:03:09 logger: small-molecule smile map:
       sm                                              smile
0    5743  C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@@]4([C@...
1   84093  CC1=CN(C(=O)NC1=O)C2CC(C(O2)COP(=O)(O)OC3CC(OC...
2  148124  CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@...
3    5757  C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3...
    prob_up       pred_type
0  0.050344  Downregulation
1  0.984405    Upregulation
2  0.998275    Upregulation
3  1.000000    Upregulation
```
:::
::::


## CLI

DeepsmirUD can also be used in shell. To know how to use, please type

```{code} shell
deepsmirud -h
```

It shows the usage of different parameters.

```{code} shell
-m, --method,
    A deep learning method. It can be any below.
    AlexNet | BiRNN | RNN | Seq2Seq |
    CNN | ConvMixer64 | DSConv | LSTMCNN |
    MobileNet | ResNet18 | ResNet50 | SEResNet
-br, --br_fpn, binary relations between small molecules and mirnas
-sm, --smile_fpn, map between small molecule IDs and their smile strings
-mir, --fasta_fp, miRNA fasta file paths
-mf, --model_fp, a model path
-o, --sv_fpn, outputting deepsmirud predictions
```

You can run it using the following code.

::::{tab-set}
:::{tab-item} Command
:sync: tab1
```{code} shell
deepsmirud -m LSTMCNN -br ./data/deepsmirud/example_data/br_sm_mirna.txt -sm ./data/deepsmirud/example_data/br_smile.txt -mir ./data/deepsmirud/example_data/ -mf ./data/deepsmirud/model/lstmcnn -o ./data/deepsmirud/out.deepsmirud
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

05/04/2025 12:08:12 logger: =>Prediction starts...
05/04/2025 12:08:12 logger: small-molecule smile map:
       sm                                              smile
0    5743  C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@@]4([C@...
1   84093  CC1=CN(C(=O)NC1=O)C2CC(C(O2)COP(=O)(O)OC3CC(OC...
2  148124  CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@...
3    5757  C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3...
[12:08:12] DEPRECATION WARNING: please use MorganGenerator
[12:08:12] DEPRECATION WARNING: please use MorganGenerator
[12:08:12] DEPRECATION WARNING: please use MorganGenerator
[12:08:12] DEPRECATION WARNING: please use MorganGenerator
    prob_up       pred_type
0  0.050344  Downregulation
1  0.984405    Upregulation
2  0.998275    Upregulation
3  1.000000    Upregulation
```
:::
::::