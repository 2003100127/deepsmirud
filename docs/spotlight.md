# Spotlight

## Sun's series work in drug discovery

Jianfeng Sun spearheaded a research plan dedicated to the AI-based discovery of small molecule therapeutics targeting non-coding RNAs and proteins, working in collaboration with both experimentalists and computational scientists across the globe. He has released 5 studies as in [Table 1](#tbl:jsun-sysbiol-work). The development of DeepdlncUD was one of the series studies of this project.

:::{table} Sun's work in drug discovery. ➵ stands for the current work.
:label: tbl:jsun-sysbiol-work
:align: center

<table border="1" cellspacing="0" cellpadding="6">
  <thead>
    <tr style="background-color:#d9d9d9;">
      <th>Field</th>
      <th>Molecule</th>
      <th>Tool name</th>
      <th>Function</th>
      <th>Technology</th>
      <th>Publication</th>
    </tr>
  </thead>
  <tbody>
    <!-- Systems Biology -->
    <tr>
      <td rowspan="3"><strong>Systems Biology</strong></td>
      <td rowspan="2">noncoding RNA</td>
      <td style="background: -webkit-linear-gradient(20deg, #09009f, #E743D9); -webkit-background-clip: text; -webkit-text-fill-color: transparent;"><strong>➵DeepsmirUD</strong></td>
      <td><em>drug discovery</em></td>
      <td>Artificial intelligence</td>
      <td>
        <a href="https://doi.org/10.3390/ijms24031878" title="DeepsmirUD">Sun et al., 2023</a>.
        <em>International Journal of Molecular Sciences</em>
      </td>
    </tr>
    <tr>
      <td><strong>DeepdlncUD</strong></td>
      <td><em>drug discovery</em></td>
      <td>Artificial intelligence</td>
      <td>
        <a href="https://doi.org/10.1016/j.compbiomed.2023.107226" title="DeepdlncUD">Sun et al., 2023</a>.
        <em>Computers in Biology and Medicine</em>
      </td>
    </tr>
    <tr>
      <td>protein</td>
      <td><strong>Drutai</strong></td>
      <td><em>drug discovery</em></td>
      <td>Artificial intelligence</td>
      <td>
        <a href="https://doi.org/10.1016/j.ejmech.2023.115500" title="Drutai">Sun et al., 2023</a>.
        <em>European Journal of Medicinal Chemistry</em>
      </td>
    </tr>
    <!-- Structural Biology -->
    <tr>
      <td rowspan="2"><strong>Structural Biology</strong></td>
      <td rowspan="2">protein</td>
      <td><strong>DeepHelicon</strong></td>
      <td><em>structural prediction</em></td>
      <td>Artificial intelligence</td>
      <td>
        <a href="https://doi.org/10.1016/j.jsb.2020.107574" title="DeepHelicon">Sun and Frishman, 2020</a>.
        <em>Journal of Structural Biology</em>
      </td>
    </tr>
    <tr>
      <td><strong>DeepTMInter</strong></td>
      <td><em>protein-protein interaction prediction</em></td>
      <td>Artificial intelligence</td>
      <td>
        <a href="https://doi.org/10.1016/j.csbj.2021.03.005" title="DeepTMInter">Sun and Frishman, 2021</a>.
        <em>Computational and Structural Biotechnology Journal</em>
      </td>
    </tr>
  </tbody>
</table>
:::

## Feature of the computational method

DeepsmirUD introduces a novel computational approach for predicting the regulatory effects of small molecules on miRNA expression by uniquely leveraging connectivity scores (<https://doi.org/10.1093/bib/bbw112>), a concept traditionally used in transcriptomic drug discovery, but here applied in an end-to-end deep learning context. Unlike prior tools that rely heavily on differential expression matrices or biological assays, DeepsmirUD infers upregulation or downregulation of miRNAs directly from molecular descriptors and sequence features, guided by connectivity-based matching between small molecules and miRNA expression profiles.

> By bypassing the need for expression count data and instead using connectivity scores as a mechanistic backbone, DeepsmirUD opens a new avenue for scalable, sequence-driven discovery of miRNA-targeting small molecules—bridging cheminformatics and systems pharmacology in a way not achieved by previous methods.


## Support for running multiple cases

In this updated version `0.1.2` in PyPI, we tuned **DeepsmirUD** to make it available to run multiple instances of predicting SM-miRNA regulation types. This should be seen as a major update because it is very important for researchers to screen large-scale regulation types with reduced oprations from their back ends.

## Runtime

**DeepsmirUD** was developed based on molecular sequences alone. It runs in a very fast speed than those supported with complex data structures and settings, making it an ideal tool for biochemical researchers.

:::{caution} Comparison
The previous version of **DeepsmirUD** runs on one regulation type prediction each time, and thus deep learning libraries will be mounted each time for multiple instances, making it time-consuming.
:::