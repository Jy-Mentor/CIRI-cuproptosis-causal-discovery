# 设置CRAN镜像
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

# 安装必要的包
if(!require(rcdk)) install.packages("rcdk")
if(!require(rJava)) install.packages("rJava")

library(rcdk)
library(rJava)

# 读取SDF文件
input_file <- "output/new_receptors/Conformer3D_COMPOUND_CID_5281515.sdf"
output_file <- "output/new_receptors/Conformer3D_COMPOUND_CID_5281515_no_h.sdf"

cat("Reading SDF file:", input_file, "\n")

# 读取分子
mols <- load.molecules(input_file)

if(length(mols) > 0) {
  mol <- mols[[1]]
  
  # 移除氢原子
  mol_no_h <- remove.hydrogens(mol)
  
  # 保存为SDF（不含氢）
  write.molecules(list(mol_no_h), output_file)
  
  cat("已生成无氢文件：", output_file, "\n")
  cat("原子数：", get.atom.count(mol_no_h), "\n")
} else {
  cat("读取失败\n")
}
