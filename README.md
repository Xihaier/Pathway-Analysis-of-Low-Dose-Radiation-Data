# pathway_analysis
Pathway analysis of low dose radiation data

<p><img src="img/method.png" title="Method" width="900"><p>

## Associated Papers
- Comprehensive analysis of gene expression profiles to radiation exposure reveals molecular signatures of low-dose radiation response [ [ArXiV](https://arxiv.org/abs/2301.01769) ] [ [2022 IEEE International Conference on Bioinformatics and Biomedicine (BIBM)](https://ieeexplore.ieee.org/abstract/document/9995607) ]
- Best poster award at 2022 Summer Symposium Adjourned Precision Medicine Applications in Radiation Oncology [ [CI4CC]([https://arxiv.org/abs/2202.11244](https://www.ci4cc.org/2022-InPerson-society-symposium)) ]

[Xihaier Luo](https://xihaier.github.io/), Balasubramanya T Nadiga, Yihui Ren, Ji Hwan Park, Wei Xu, Shinjae Yoo


## Dependencies
- python 3
- PyTorch
- matplotlib


## Installation

- Install PyTorch and other dependencies

- Clone this repo:

```
git clone https://github.com/Xihaier/Near-Term-Climate-Prediction-BDL
```


## Dataset

The dataset used have been uploaded to Google Drive and can be downloaded with corresponding links.

Link: https://zenodo.org/record/6822275#.YzMpJOyZMeY


## Model
<p><img src="img/problem.png" title="problem" width="500"><p>
  
We ran extensive models and shared some of the best ones here, which can be divided into two categories: deterministic models and Bayesian models. For example, if a deterministic model, such as DenseNet, is chosen, one should run 

```bash
cd src
python main_DL.py
```

If you are interested in ConvLSTM and its variants, we have included our implementation in the appendix folder.

```bash
cd appendix
python main.py
```

## Citation

If you find this repo useful for your research, please consider to cite:

```latex
@article{luo2022bayesian,
  title={A Bayesian Deep Learning Approach to Near-Term Climate Prediction},
  author={Luo, Xihaier and Nadiga, Balasubramanya T and Ren, Yihui and Park, Ji Hwan and Xu, Wei and Yoo, Shinjae},
  journal={arXiv preprint arXiv:2202.11244},
  year={2022}
}
```

## Questions

For any questions or comments regarding this paper, please contact Xihaier Luo via [xluo@bnl.gov](mailto:xluo@bnl.gov).

