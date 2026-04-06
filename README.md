# multiPrimerEvaluation
tNGS 单管超多重引物评估工具 - 用于评估多重引物在同一管反应体系中的相互作用风险，rust语言实现;
> **划重点**：<br>
> 支持多重一管组合评估，上千重同样支持；<br>
> tNGS 超多重组合逻辑：作减法，迭代优化;<br>
> badness值评估核心逻辑为 **确保引物3'端的结合活性**<br>
> 引物独立设计，推荐使用**RUCS设计物种保守特异引物**<br>

## 使用方法

```bash
git clone https://github.com/yi1873/multiPrimerEvaluation.git
cd multiPrimerEvaluation
cargo build --release

./target/release/multiPrimerEvaluation -h
## v0.2.0版本 支持多线程并行处理，默认4线程

## input: primer.tsv
./target/release/multiPrimerEvaluation -i test/primer.tsv -o test/primer.badnessDisk.inOneTube.tsv
## input: primer.fa
./target/release/multiPrimerEvaluation -i test/primer.fa  -o test/primer.badnessDisk.inOneTube.tsv 
```

## 项目简介

本项目用于评估多个引物在同一反应管中可能存在的相互作用风险，通过计算引物间的badness值并进行风险等级评估，帮助避免引物二聚体、非特异性扩增等问题。

**注意**：本项目现在提供Rust实现版本，性能更高，功能更全面。

## 核心功能

- **反向互补检测**：识别引物间形成的反向互补序列
- **Badness计算**：基于反向互补序列长度、距离等参数量化引物相互作用风险
- **风险等级评估**：将badness值转换为低/中/高风险三个等级
- **错配容忍**：支持指定最大错配数（默认2个）
- **阈值过滤**：只输出badness值大于阈值的引物对（默认10）

## Badness计算方法

### 计算公式

```
badness(Pa,Pb) = P × (2^len) / ((d1+1) × (d2+1))
```

**参数说明：**
- `len`：引物间反向互补子序列的长度
- `d1`：反向互补序列到Primer_1末端的距离
- `d2`：反向互补序列到Primer_2末端的距离
- `P`：权重因子
  - d1和d2均为0时，P = 100（提升2个数量级）
  - d1和d2中只有一个为0时，P = 10（提升1个数量级）
  - d1和d2都不为0时，P = 1

### 风险等级划分

- **高风险**：log10(badness) ≥ 3.5（引物末端形成5bp以上反向互补）
- **中风险**：3 ≤ log10(badness) < 3.5
- **低风险**：log10(badness) < 3

### 参考文献 【思路来源，实现逻辑有加权处理】

- Xie NG, Wang MX, Song P, Mao S, Wang Y, Yang Y, Luo J, Ren S, Zhang DY. Designing highly multiplex PCR primer sets with Simulated Annealing Design using Dimer Likelihood Estimation (SADDLE). Nat Commun. 2022 Apr 11;13(1):1881. doi: 10.1038/s41467-022-29500-4. PMID: 35410464; PMCID: PMC9001684.

## Citation 📚

If you use this tool in your research, please cite:

```bibtex
@software{multiPrimerEvaluation,
  title = {multiPrimerEvaluation: High-performance primer badness evaluation in Rust},
  author = {Liang Xiangzhi},
  year = {2026},
  url = {https://github.com/yi1873/multiPrimerEvaluation}
}
```


