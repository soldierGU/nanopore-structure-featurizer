"""
npstructfeat

Nanopore Structure Featurizer 的核心 Python 包。

这个包主要用于：
1. 读取 PDB/mmCIF 结构文件；
2. 提取纳米孔蛋白的残基级结构特征；
3. 筛选候选孔道内壁残基；
4. 聚合得到纳米孔整体结构特征表。

注意：
scripts/ 目录只负责运行流程；
真正可复用的函数应放在 src/npstructfeat/ 下面。
"""


__version__ = "0.1.0"