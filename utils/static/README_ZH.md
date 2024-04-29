# IAGS_AUTO
基于IAGS的自动化流程工具

[English Version](../../README.md)

## 下载

可以通过conda，或者下载源码来使用此工具。

### conda下载（推荐）

+ 创建虚拟环境
   ```shell
   conda create -n iags_auto python=3.9
   ```
+ 下载
   ```shell
   conda install -c gurobi -c conda-forge -c huntguo iags_auto
   ```

### 源码使用

+ 下载源代码
   ```shell
   wget https://codeload.github.com/99gloom/IAGS_AUTO/zip/refs/heads/main
   ```
+ 下载mono
   ```shell
   sudo apt install mono-devel
   ```

*注意：两种方法下载后，都需要激活[gurobi](https://www.gurobi.com)，这里我们提供了一份[帮助文档](gurobi_ZH.md)，帮助您获取许可*。

## 使用

### 文件

IAGS_AUTO工具需要三种文件：各物种gff文件、orthogroup.tsv文件，以及species.tree文件，并将这三种文件放置于同一文件夹下。
其中gff和orthogroup.tsv文件与之前的脚本[processDrimm](https://github.com/99gloom/processDrimm) 所需的相同。

下面将逐一介绍：

1. gff文件：格式与[MCScanX](https://github.com/wyp1125/MCScanx) 输入的相同。文件内只有四列，分别是染色体名称、基因名、基因起始坐标、基因结束坐标。格式如下：  
    ```text
    sp_name  gene_name  starting_position  ending_position
    ```

2. Orthogroups.tsv：[OrthoFinder](https://github.com/davidemms/OrthoFinder) 的输出文件

3. species.tree：WGD-Newick格式。本质为Newick格式的修改版本，在树中的WGD位点，添加了`[WGD]`标志。下图中红色点位为WGD标志。（结尾有无`';'`号都可以）  
  ![图中红色点位为WGD标志](img.png)

### 用法

<table>
<tr>
<th width="120">命令</th>
<th width="200">参数</th>
<th width="400">说明</th>
</tr>
<tr>
<td>-f, --filepath</td>
<td>./file_dir</td>
<td>三种需要的文件存放的目录</td>
</tr>
<tr>
<td>-c, --cycleLength</td>
<td>默认为20</td>
<td>控制共线块的连续性</td>
</tr>
<tr>
<td>-d, --dustLength</td>
<td>默认为所有物种拷贝数+1</td>
<td>控制基因家族大小的上限，会将同源基因个数中超过这个上限的基因家族进行过滤</td>
</tr>
<tr>
<td rowspan="2">-s, --shape</td>
<td>"s"（默认）</td>
<td rowspan="2">染色体形状，s为线性染色体，c为环状染色体</td>
</tr>
<tr>
<td>"c"</td>
</tr>
<tr>
<td rowspan="2">-m, --model</td>
<td>"manual"</td>
<td rowspan="2">默认为None。当用户需要手动指定outgroup时，首先使用"manual"模式生成节点计算顺序文件 "model_and_outgroup.txt"，修改其使用的outgroup信息。再使用"continue"模式，基于上一步修改的问题件生成结果。</td>
</tr>
<tr>
<td>"continue"</td>
</tr>
<tr>
<td>"--dotplot"</td>
<td>-</td>
<td>生成各物种两两之间的dotplot图</td>
</tr>
<tr>
<td>"--expand"</td>
<td>-</td>
<td>通过图扩充算法增加共线块覆盖度</td>
</tr>
<tr>
<td rowspan="2">--check</td>
<td>"yes" (Default)</td>
<td rowspan="2">按拷贝数比例过滤共线块后，空染色体占比大于30%时是否停止程序（共线块质量低时停止程序）</td>
</tr>
<tr>
<td>"no"</td>
</tr>
</table>

### 结果

运行后会在运行目录下生成Result文件夹，内部主要有子文件夹，按生成的顺序分别为Tree_File，Process_Drimm，IAGS。
用户主要关注的文件为IAGS以及Tree_File中的model_and_outgroup.txt。IAGS为最终生成的祖先基因组结果，以及染色体图。
model_and_outgroup.txt为各祖先节点所用外族信息，如需手动指定则需要额外处理。

下面将详细介绍各个文件的作用：
+ Tree_File
  + species.ratio和all.ratio：所有物种的拷贝数信息；
  + Evolutionary_tree.txt：进化树树形以及所有节点分布；
  + model_and_outgroup.txt：各个祖先节点计算信息，格式为：当前计算的祖先节点:IAGS计算模型:孩子节点:外族。如果模型为GMP或MultiCopyGMP，则孩子节点有两个。尤其是MultiCopyGMP模型中，如果其外族的拷贝数不足，则会对外族染色体手动加倍以计算，以"*N"表示。
+ Process_Drimm  
  本质为[processDrimm](https://github.com/99gloom/processDrimm) 的自动化过程。
  + Process_OrthoFind：通过物种GFF和orthogroup.tsv文件对基因编码，生成基因序列；
  + Drimm_Synteny_Output：运行DRIMM后的原始结果；
  + Drimm_Blocks：对DRIMM运行后的原始结果进行LCS，以便下游分析；
  + Final_Blocks：将Drimm_Blocks按比例过滤，生成可供IAGS运行的block。
+ IAGS
  + 各祖先节点名：祖先节点详细信息，包括计算的block，CRB ratio评估等；
  + painting：祖先基因组染色体图，其中Painting_start_point.txt记录了画图的基点祖先；
  + shufflingEvents.txt：物种fission和fusion信息。
  
### 示例

+ 1.快速开始
  ```shell
  iags_auto -f ./example
  ```

+ 2.指定参数
  ```shell
  iags_auto -f ./example -c 60 -d 12 -s s
  ```
+ 3.通过图扩充算法运行
  
  ```shell
  python IAGS_AUTO -f ./example --expand
  ```
  
+ 4.手动指定外族
  ```shell
  iags_auto -f ./example -m manual
  ```
  随后更改 Result/Tree_File/model_and_outgroup.txt中的外族后继续运行。
  ```shell
  iags_auto -f ./example -m continue
  ```
  
+ 5.绘制dotplot
  ```shell
  iags_auto -f ./example --dotplot
  ```
  使用此命令，会在Result中额外生成一个Dotplot文件夹，将结果存放其中。

+ 6.不进行染色体质量检验
  ```shell
  iags_auto -f ./example --check no
  ```

如果您选择直接下载源码使用，则需要将上述命令中`iags_auto`替换为`python IAGS_ATUO.py`，通过直接调用python文件启动。

## 其他

由于IAGS基于[gurobi](https://www.gurobi.com)进行整数优化，所以此工具需要用户自行下载并激活gurobi的license，这里我们提供了的[帮助文档](gurobi_ZH.md)帮助用户安装激活gurobi。












 