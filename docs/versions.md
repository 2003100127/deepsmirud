# Versions

## Change

::::{grid} 2 2 2 2

| Version |                                   Date                                   |
|:-------:|:------------------------------------------------------------------------:|
| `0.1.1` |   ![](https://img.shields.io/badge/past_released-March._2023-red.svg)    |
| `0.1.2` | ![](https://img.shields.io/badge/latest_released-April._2025-green.svg)  |

::::

## Before `0.1.2`

The version `0.1.1` of **DeepsmirUD** was used to infer a small molecule-miRNA regulation type per instance only. This process is simplified as follows.

### Small molecule 

> Estradiol (CID: 5757)

```shell
C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O
```

### Target

> MIMAT0000066

```
>MIMAT0000066 (hsa-let-7e-5p)
UGAGGUAGGAGGUUGUAUAGUU
```

### Inference

```shell
deepsmirud -m LSTMCNN -sm deepsmirud/data/example/5757.txt -mir deepsmirud/data/example/MIMAT0000066.fasta -mf /the/path/you/prefer/model/lstmcnn -o ./out.deepsmirud
```


## Tensorflow

::::{caution} Major change in deep learing libraries
In Keras 3, the way to call models becomes completely different from that in Keras 2.

:::{table} The way to use models.
:label: tbl1
:align: center

| Item              | TF Keras 2                 | TF Keras 3                       |
|-------------------|----------------------------|----------------------------------|
| model format      | `.h5` or TF `SavedModel`   | `.keras` (new) + `.h5`  (legacy) |
| `subclass` Models | supported via `SavedModel` | `TFSMLayer` for loading          |

:::

::::

### Comparison

Our experience suggests that changes in Python versions have a smaller impact compared to changes in TensorFlow versions.

:::{table} The way to use models.
:label: tbl2
:align: center

| Python  | Tensorflow verisons        | Supported versions     |
|---------|----------------------------|------------------------|
| `3.11`  | 🌟`2.14`🌟, `2.15`, `2.16` | TF Keras 2, TF Keras 3 |
| `3.12`  | `2.17`, `2.18`, `2.19`     | only TF Keras 3        |

:::